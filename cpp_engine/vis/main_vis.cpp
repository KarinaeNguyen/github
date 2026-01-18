// main_vis.cpp (based on user's working version; adds square ceiling rail + long-run plot fixes)
// - Keeps fire_center locked to simulation truth via refresh_obs() (Observation.hotspot_pos_m_*)
// - Does NOT require any new Observation fields (no nozzle truth) or new Simulation methods
// - Renders a model-backed square ceiling rail (CeilingRail) that auto-resizes with rack_half
//   and is fixed 30 cm below ceiling (configurable via UI)
// - Fixes "plots look broken" by:
//     * driving simTime from sim.time_s() (authoritative sim clock)
//     * showing Samples + time range
//     * forcing X axis limits to the current window [t0, t1] for each plot

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstdint>

#include "Simulation.h"

// Step 1: model-backed ceiling rail (no UI dependencies)
#include "../world/ceiling_rail.h"
#include "../world/rail_mounted_nozzle.h"

#include "imgui.h"
// ---- Docking compatibility shim (older ImGui builds do not define docking flags/APIs)
#ifndef ImGuiConfigFlags_DockingEnable
#define VFEP_NO_IMGUI_DOCKING 1
#endif
#include "implot.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "imgui_internal.h"  // DockSpaceOverViewport / DockBuilder

// Platform GL headers: on Windows, <GL/gl.h> requires Windows types/macros (APIENTRY/WINGDIAPI).
// Include <windows.h> first to avoid syntax errors in the Windows SDK gl.h.
#ifdef _WIN32
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <windows.h>
#endif

#include <GLFW/glfw3.h>

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif

// ============================================================
// Logo Loading & Management (Simple PPM format for zero dependencies)
// ============================================================

struct LogoTexture {
    GLuint texture_id = 0;
    int width = 0;
    int height = 0;
    bool loaded = false;
};

// Load JPEG logo (we'll use a simple fallback approach)
// For now, we'll create a procedural logo if file not found
static LogoTexture load_logo_texture(const char* logo_path) {
    LogoTexture logo;
    
    // Try to load from file using simple approach
    FILE* f = fopen(logo_path, "rb");
    if (!f) {
        // Fallback: create simple procedural logo (green rectangle with text)
        // We'll render it via ImGui instead of texture
        logo.loaded = false;
        return logo;
    }
    
    // For JPG support, you can add stb_image.h later
    // For now, just close and return false
    fclose(f);
    logo.loaded = false;
    return logo;
}

static void glfw_error_callback(int error, const char* description) {
    std::fprintf(stderr, "GLFW Error %d: %s\n", error, description ? description : "(null)");
}

static int fail(const char* msg) {
    std::fprintf(stderr, "FATAL: %s\n", msg ? msg : "(null)");
    std::fprintf(stderr, "\n");
    return EXIT_FAILURE;
}

// ============================================================
// Minimal Phase-1 3D Twin (fixed-pipeline, deterministic, no assets)
// Adds: HUD overlay + visualization toggles + docking
// ============================================================

struct Vec3f { float x, y, z; };

static Vec3f v3(float x, float y, float z) { return {x,y,z}; }

// Conversions between vis-local float vectors and model double vectors.
static vfep::world::Vec3d to_v3d(Vec3f v) {
    return vfep::world::Vec3d{ (double)v.x, (double)v.y, (double)v.z };
}

static Vec3f to_v3f(const vfep::world::Vec3d& v) {
    return v3((float)v.x, (float)v.y, (float)v.z);
}

static Vec3f add(Vec3f a, Vec3f b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
static Vec3f sub(Vec3f a, Vec3f b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
static Vec3f mul(Vec3f a, float s)  { return {a.x*s, a.y*s, a.z*s}; }

static float clampf(float x, float lo, float hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

static float dot(Vec3f a, Vec3f b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static Vec3f cross(Vec3f a, Vec3f b) { return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x }; }
static float len(Vec3f a) { return std::sqrt(dot(a,a)); }
static Vec3f norm(Vec3f a) {
    float l = len(a);
    return (l > 1e-6f) ? mul(a, 1.0f/l) : v3(0,0,0);
}

// Step 2 staging: axis-angle rotation (Rodrigues)
static Vec3f rotate_axis_angle(Vec3f v, Vec3f axis_unit, float ang_rad) {
    Vec3f a = norm(axis_unit);
    if (len(a) < 1e-6f) return v;
    const float c = std::cos(ang_rad);
    const float s = std::sin(ang_rad);
    // Rodrigues: v' = v*c + (a×v)*s + a*(a·v)*(1-c)
    return add(
        add(mul(v, c), mul(cross(a, v), s)),
        mul(a, dot(a, v) * (1.0f - c))
    );
}

static void set_perspective(float fovy_deg, float aspect, float znear, float zfar) {
    // OpenGL fixed pipeline expects column-major matrix.
    const float fovy_rad = fovy_deg * 3.1415926535f / 180.0f;
    const float f = 1.0f / std::tan(0.5f * fovy_rad);

    float m[16] = {};
    m[0]  = f / aspect;
    m[5]  = f;
    m[10] = (zfar + znear) / (znear - zfar);
    m[11] = -1.0f;
    m[14] = (2.0f * zfar * znear) / (znear - zfar);

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(m);
}

static void look_at(Vec3f eye, Vec3f center, Vec3f up) {
    // Minimal lookAt for fixed pipeline.
    Vec3f fwd = sub(center, eye);
    float fl = len(fwd);
    if (fl > 1e-6f) fwd = mul(fwd, 1.0f / fl);

    float ul = len(up);
    if (ul > 1e-6f) up = mul(up, 1.0f / ul);

    Vec3f s = cross(fwd, up);
    float sl = len(s);
    if (sl > 1e-6f) s = mul(s, 1.0f / sl);

    Vec3f u = cross(s, fwd);

    float m[16] = {
        s.x,  u.x,  -fwd.x, 0.0f,
        s.y,  u.y,  -fwd.y, 0.0f,
        s.z,  u.z,  -fwd.z, 0.0f,
        0.0f, 0.0f, 0.0f,   1.0f
    };

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(m);
    glTranslatef(-eye.x, -eye.y, -eye.z);
}

static void draw_wire_box(Vec3f c, Vec3f half) {
    const float x0 = c.x - half.x, x1 = c.x + half.x;
    const float y0 = c.y - half.y, y1 = c.y + half.y;
    const float z0 = c.z - half.z, z1 = c.z + half.z;

    glBegin(GL_LINES);
    // bottom
    glVertex3f(x0,y0,z0); glVertex3f(x1,y0,z0);
    glVertex3f(x1,y0,z0); glVertex3f(x1,y0,z1);
    glVertex3f(x1,y0,z1); glVertex3f(x0,y0,z1);
    glVertex3f(x0,y0,z1); glVertex3f(x0,y0,z0);
    // top
    glVertex3f(x0,y1,z0); glVertex3f(x1,y1,z0);
    glVertex3f(x1,y1,z0); glVertex3f(x1,y1,z1);
    glVertex3f(x1,y1,z1); glVertex3f(x0,y1,z1);
    glVertex3f(x0,y1,z1); glVertex3f(x0,y1,z0);
    // verticals
    glVertex3f(x0,y0,z0); glVertex3f(x0,y1,z0);
    glVertex3f(x1,y0,z0); glVertex3f(x1,y1,z0);
    glVertex3f(x1,y0,z1); glVertex3f(x1,y1,z1);
    glVertex3f(x0,y0,z1); glVertex3f(x0,y1,z1);
    glEnd();
}

static void draw_solid_box(Vec3f c, Vec3f half) {
    const float x0 = c.x - half.x, x1 = c.x + half.x;
    const float y0 = c.y - half.y, y1 = c.y + half.y;
    const float z0 = c.z - half.z, z1 = c.z + half.z;

    glBegin(GL_QUADS);
    // +Z
    glVertex3f(x0,y0,z1); glVertex3f(x1,y0,z1); glVertex3f(x1,y1,z1); glVertex3f(x0,y1,z1);
    // -Z
    glVertex3f(x1,y0,z0); glVertex3f(x0,y0,z0); glVertex3f(x0,y1,z0); glVertex3f(x1,y1,z0);
    // +X
    glVertex3f(x1,y0,z1); glVertex3f(x1,y0,z0); glVertex3f(x1,y1,z0); glVertex3f(x1,y1,z1);
    // -X
    glVertex3f(x0,y0,z0); glVertex3f(x0,y0,z1); glVertex3f(x0,y1,z1); glVertex3f(x0,y1,z0);
    // +Y
    glVertex3f(x0,y1,z1); glVertex3f(x1,y1,z1); glVertex3f(x1,y1,z0); glVertex3f(x0,y1,z0);
    // -Y
    glVertex3f(x0,y0,z0); glVertex3f(x1,y0,z0); glVertex3f(x1,y0,z1); glVertex3f(x0,y0,z1);
    glEnd();
}

static void draw_line(Vec3f a, Vec3f b) {
    glBegin(GL_LINES);
    glVertex3f(a.x,a.y,a.z);
    glVertex3f(b.x,b.y,b.z);
    glEnd();
}

static void temp_to_color(float tempC, float& r, float& g, float& b) {
    // Conservative ramp. 24C neutral -> hotter -> red.
    float t = clampf((tempC - 24.0f) / (120.0f - 24.0f), 0.0f, 1.0f);
    // gray -> yellow -> orange -> red
    r = 0.25f + 0.75f * t;
    g = 0.25f + 0.50f * (1.0f - std::abs(2.0f*t - 1.0f)); // peak mid
    b = 0.25f * (1.0f - t);
}

static const char* suppression_regime_text(int r) {
    switch (r) {
        case 0: return "None";
        case 1: return "Ineffective";
        case 2: return "Marginal";
        case 3: return "Effective";
        case 4: return "Overkill";
        default: return "Unknown";
    }
}

static float fire_scale_from_HRR_W(double hrrW) {
    // cube-root scaling for stability (use kW ref)
    const float hrr_kW = (float)(hrrW * 0.001);
    const float ref_kW = 1000.0f; // 1 MW reference
    float s = std::pow(std::max(hrr_kW, 0.0f) / ref_kW, 1.0f / 3.0f);
    return clampf(s, 0.10f, 2.00f);
}

// Draw a simple cone aligned to dir_unit in world space (apex at nozzle).
static void draw_cone_world(Vec3f apex, Vec3f dir_unit, float length_m, float radius_m, int slices = 16) {
    Vec3f d = norm(dir_unit);
    if (len(d) < 1e-6f || length_m <= 1e-4f || radius_m <= 1e-4f) return;

    // Orthonormal basis with z = d
    Vec3f up = (std::abs(d.y) < 0.9f) ? v3(0,1,0) : v3(1,0,0);
    Vec3f x = norm(cross(up, d));
    Vec3f y = cross(d, x);

    Vec3f base_center = add(apex, mul(d, length_m));

    glBegin(GL_TRIANGLES);
    for (int i = 0; i < slices; ++i) {
        const float a0 = (2.0f * 3.1415926535f * (float)i) / (float)slices;
        const float a1 = (2.0f * 3.1415926535f * (float)(i+1)) / (float)slices;

        Vec3f p0 = add(base_center, add(mul(x, radius_m * std::cos(a0)), mul(y, radius_m * std::sin(a0))));
        Vec3f p1 = add(base_center, add(mul(x, radius_m * std::cos(a1)), mul(y, radius_m * std::sin(a1))));

        glVertex3f(apex.x, apex.y, apex.z);
        glVertex3f(p0.x,   p0.y,   p0.z);
        glVertex3f(p1.x,   p1.y,   p1.z);
    }
    glEnd();
}

// Very simple arrow: shaft + line head.
static void draw_arrow(Vec3f origin, Vec3f dir_unit, float length_m) {
    Vec3f d = norm(dir_unit);
    if (len(d) < 1e-6f || length_m <= 1e-4f) return;

    Vec3f tip = add(origin, mul(d, length_m));

    Vec3f up = (std::abs(d.y) < 0.9f) ? v3(0,1,0) : v3(1,0,0);
    Vec3f x = norm(cross(up, d));
    Vec3f y = cross(d, x);

    const float head_len = length_m * 0.18f;
    const float head_w   = length_m * 0.06f;

    Vec3f h0 = add(tip, add(mul(d, -head_len), mul(x,  head_w)));
    Vec3f h1 = add(tip, add(mul(d, -head_len), mul(x, -head_w)));
    Vec3f h2 = add(tip, add(mul(d, -head_len), mul(y,  head_w)));
    Vec3f h3 = add(tip, add(mul(d, -head_len), mul(y, -head_w)));

    draw_line(origin, tip);
    draw_line(tip, h0);
    draw_line(tip, h1);
    draw_line(tip, h2);
    draw_line(tip, h3);
}

// Ray vs axis-aligned box intersection (slab method).
// Returns true if intersects; t_hit is distance along ray to first hit (>=0).
static bool ray_aabb_intersect(Vec3f ro, Vec3f rd_unit, Vec3f box_center, Vec3f box_half, float& t_hit) {
    Vec3f rd = rd_unit;
    // Avoid divide-by-zero; treat near-zero components as huge inv.
    auto inv = [&](float v) -> float { return (std::abs(v) > 1e-8f) ? (1.0f / v) : 1e30f; };

    float tmin = -1e30f;
    float tmax =  1e30f;

    const float bmin_x = box_center.x - box_half.x;
    const float bmax_x = box_center.x + box_half.x;
    const float bmin_y = box_center.y - box_half.y;
    const float bmax_y = box_center.y + box_half.y;
    const float bmin_z = box_center.z - box_half.z;
    const float bmax_z = box_center.z + box_half.z;

    float tx1 = (bmin_x - ro.x) * inv(rd.x);
    float tx2 = (bmax_x - ro.x) * inv(rd.x);
    tmin = std::max(tmin, std::min(tx1, tx2));
    tmax = std::min(tmax, std::max(tx1, tx2));

    float ty1 = (bmin_y - ro.y) * inv(rd.y);
    float ty2 = (bmax_y - ro.y) * inv(rd.y);
    tmin = std::max(tmin, std::min(ty1, ty2));
    tmax = std::min(tmax, std::max(ty1, ty2));

    float tz1 = (bmin_z - ro.z) * inv(rd.z);
    float tz2 = (bmax_z - ro.z) * inv(rd.z);
    tmin = std::max(tmin, std::min(tz1, tz2));
    tmax = std::min(tmax, std::max(tz1, tz2));

    if (tmax < 0.0f) return false;       // box behind ray
    if (tmin > tmax) return false;

    t_hit = (tmin >= 0.0f) ? tmin : tmax; // if inside box, take exiting hit
    return t_hit >= 0.0f;
}

struct VisualUIState {
    bool show_hud = true;
    bool show_controls = true;
    bool show_plots = true;

    bool draw_warehouse = true;
    bool draw_rack = true;
    bool draw_fire = true;
    bool draw_fire_sectors = true;
    bool draw_draft = true;
    bool draw_nozzle = true;
    bool draw_spray = true;
    bool draw_hit_marker = true;

    bool draw_ceiling_rail = true;
};

static void plot_line_with_xlimits(const char* title,
                                  const char* label,
                                  const double* xs,
                                  const double* ys,
                                  int count,
                                  double t0,
                                  double t1)
{
    if (count <= 1)
        return;

    if (ImPlot::BeginPlot(title)) {

        // --- X-axis handling (robust across ImPlot versions) ---
#if defined(ImAxis_X1)
        // ImPlot >= 0.16
        ImPlot::SetupAxisLimits(ImAxis_X1, t0, t1, ImGuiCond_Always);
#elif defined(ImPlotAxis_X1)
        // Transitional versions
        ImPlot::SetupAxisLimits(ImPlotAxis_X1, t0, t1, ImGuiCond_Always);
#else
        // Very old ImPlot: DO NOT set limits (auto-fit fallback)
        // This avoids calling deprecated/nonexistent APIs.
#endif

        // --- Plot data ---
        ImPlot::PlotLine(label, xs, ys, count);

        ImPlot::EndPlot();
    }
}

int main(int argc, char** argv) {
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit()) return fail("glfwInit failed");

    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

    GLFWwindow* window = glfwCreateWindow(1280, 720, "VFEP Visualizer", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return fail("glfwCreateWindow failed");
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // vsync

    // Validate OpenGL context exists.
    const GLubyte* gl_version = glGetString(GL_VERSION);
    if (!gl_version) {
        glfwDestroyWindow(window);
        glfwTerminate();
        return fail("OpenGL context validation failed (glGetString(GL_VERSION) returned null)");
    }
    std::fprintf(stderr, "OpenGL Vendor:   %s\n", glGetString(GL_VENDOR));
    std::fprintf(stderr, "OpenGL Renderer: %s\n", glGetString(GL_RENDERER));
    std::fprintf(stderr, "OpenGL Version:  %s\n", gl_version);

    bool imgui_ctx = false;
    bool implot_ctx = false;
    bool imgui_glfw = false;
    bool imgui_gl3 = false;

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    imgui_ctx = true;

    ImGuiIO& io = ImGui::GetIO();
#ifndef VFEP_NO_IMGUI_DOCKING
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
#endif

    // Increase font size for better readability (19pt for better visibility)
    io.Fonts->AddFontFromFileTTF("C:/Windows/Fonts/segoeui.ttf", 19.0f);
    if (io.Fonts->Fonts.Size == 1) {
        // If font loading fails, use default with scaling
        ImGui::GetStyle().ScaleAllSizes(1.5f);
    }

    ImPlot::CreateContext();
    implot_ctx = true;

    ImGui::StyleColorsDark();

    if (!ImGui_ImplGlfw_InitForOpenGL(window, true)) {
        if (implot_ctx) ImPlot::DestroyContext();
        if (imgui_ctx) ImGui::DestroyContext();
        glfwDestroyWindow(window);
        glfwTerminate();
        return fail("ImGui_ImplGlfw_InitForOpenGL failed");
    }
    imgui_glfw = true;

    if (!ImGui_ImplOpenGL3_Init(glsl_version)) {
        if (imgui_glfw) ImGui_ImplGlfw_Shutdown();
        if (implot_ctx) ImPlot::DestroyContext();
        if (imgui_ctx) ImGui::DestroyContext();
        glfwDestroyWindow(window);
        glfwTerminate();
        return fail("ImGui_ImplOpenGL3_Init failed");
    }
    imgui_gl3 = true;

    vfep::Simulation sim;
    bool running = false;

    // --- CLI flags ---
    bool start_calib = false;
    for (int i = 1; i < argc; ++i) {
        if (argv[i] && std::string(argv[i]) == "--calib") {
            start_calib = true;
        }
    }

    VisualUIState ui;

    // --- Phase-1 3D Twin camera (deterministic, ImGui-controlled) ---
    float cam_yaw_deg   = 35.0f;
    float cam_pitch_deg = 20.0f;
    float cam_dist      = 8.0f;
    Vec3f cam_target    = v3(0.0f, 1.2f, 0.0f);

    // --- Phase-1 scene layout (meters, simple boxes) ---
    Vec3f warehouse_half = v3(6.0f, 3.0f, 6.0f);
    Vec3f rack_center    = v3(0.0f, 1.0f, 0.0f);
    Vec3f rack_half      = v3(0.6f, 1.0f, 0.4f);
    Vec3f fire_center    = v3(0.0f, 0.6f, 0.7f);

    // --- Ceiling rail params (UI -> model config) ---
    float rail_ceiling_drop_m = 0.30f; // 30 cm below ceiling
    float rail_margin_m       = 0.25f; // rail is wider than rack footprint

    // --- Ceiling rail (model-backed) ---
    vfep::world::CeilingRail       ceiling_rail;
    vfep::world::CeilingRailInputs ceiling_rail_in;
    vfep::world::CeilingRailConfig ceiling_rail_cfg;

    // --- Step 2: rail-mounted nozzle (model-backed, kinematics only) ---
    vfep::world::RailMountedNozzle       rail_nozzle;
    vfep::world::RailMountedNozzle::Config rail_nozzle_cfg;
    vfep::world::RailMountedNozzle::Inputs rail_nozzle_in;

    // UI parameter for nozzle drop below the rail
    float nozzle_drop_from_rail_m = 0.15f;

    // --- Spray/nozzle parameters ---
    Vec3f nozzle_pos     = v3(-2.0f, 1.5f, -2.0f);
    Vec3f nozzle_dir     = v3(0.7f, -0.15f, 0.7f);
    float mdot_ref       = 0.15f;
    float spray_L0       = 0.6f;
    float spray_L1       = 3.2f;
    float spray_R0       = 0.10f;
    float spray_R1       = 0.28f;
    float spray_max_len  = 8.0f;

    // --- Step 2 staging: visualization-only nozzle controls (no Simulation coupling) ---
    float viz_nozzle_s_0_1         = 0.25f;   // param along rail (future: CeilingRail)
    float viz_nozzle_pan_deg       = 0.0f;    // yaw
    float viz_nozzle_tilt_deg      = 0.0f;    // pitch
    bool  viz_override_nozzle_pose = true;


    float hit_marker_base = 0.06f;
    float hit_marker_gain = 0.20f;

    Vec3f draft_vel_mps  = v3(0.0f, 0.0f, 0.0f);
    float draft_arrow_scale = 0.7f;

    float draft_deflect_gain = 0.35f;

    double dt = 0.05;
    double simTime = 0.0;

    double wall_prev = glfwGetTime();
    double accum_s = 0.0;

    int last_substeps = 0;
    bool dropped_accum = false;

    // History buffers
    std::vector<double> t_hist, T_hist, HRR_hist, O2_hist, EffExp_hist, KD_hist, KDTarget_hist;
    t_hist.reserve(20000);
    T_hist.reserve(20000);
    HRR_hist.reserve(20000);
    O2_hist.reserve(20000);
    EffExp_hist.reserve(20000);
    KD_hist.reserve(20000);
    KDTarget_hist.reserve(20000);

    constexpr size_t kMaxHistory = 200000;
    constexpr size_t kTrimChunk  = 10000;
    constexpr int kPlotWindowN   = 5000;

    auto trim_history_if_needed = [&]() {
        if (t_hist.size() <= kMaxHistory) return;
        const size_t drop = std::min(kTrimChunk, t_hist.size());
        auto erase_front = [&](std::vector<double>& v) {
            v.erase(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(drop));
        };
        erase_front(t_hist);
        erase_front(T_hist);
        erase_front(HRR_hist);
        erase_front(O2_hist);
        erase_front(EffExp_hist);
        erase_front(KD_hist);
        erase_front(KDTarget_hist);
    };

    auto push_sample = [&](double t, const vfep::Observation& o) {
        t_hist.push_back(t);
        T_hist.push_back(o.T_K);
        HRR_hist.push_back(o.HRR_W);
        O2_hist.push_back(o.O2_volpct);
        EffExp_hist.push_back(o.effective_exposure_kg);
        KD_hist.push_back(o.knockdown_0_1);

        double kd_t = 0.0;
        for (int i = 0; i < vfep::Observation::kNumSuppressionSectors; ++i) {
            kd_t += o.sector_knockdown_target_0_1[i];
        }
        kd_t /= (double)vfep::Observation::kNumSuppressionSectors;
        KDTarget_hist.push_back(kd_t);

        trim_history_if_needed();
    };

    vfep::Observation last_obs = sim.observe();
    simTime = sim.time_s();
    push_sample(simTime, last_obs);

    // Keep UI fire marker locked to simulation truth.
    fire_center = v3((float)last_obs.hotspot_pos_m_x,
                     (float)last_obs.hotspot_pos_m_y,
                     (float)last_obs.hotspot_pos_m_z);

    // One canonical refresh point so we cannot miss an update site.
    auto refresh_obs = [&]() {
        last_obs = sim.observe();
        simTime = sim.time_s();
        fire_center = v3((float)last_obs.hotspot_pos_m_x,
                         (float)last_obs.hotspot_pos_m_y,
                         (float)last_obs.hotspot_pos_m_z);
    };

    // --- Auto-start calibration mode (deterministic) if requested ---
    if (start_calib) {
        sim.enableCalibrationMode(true);
        running = false;
        refresh_obs();
        accum_s = 0.0;

        t_hist.clear(); T_hist.clear(); HRR_hist.clear(); O2_hist.clear();
        EffExp_hist.clear(); KD_hist.clear(); KDTarget_hist.clear();

        push_sample(simTime, last_obs);
        last_substeps = 0;
        dropped_accum = false;
    }

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // --- advance sim (wall-time accumulator) ---
        const double wall_now = glfwGetTime();
        double wall_dt = wall_now - wall_prev;
        wall_prev = wall_now;

        wall_dt = std::clamp(wall_dt, 0.0, 0.1);
        dt = std::clamp(dt, 0.001, 1.0);

        bool advanced_this_frame = false;

        if (running && !sim.isConcluded()) {
            accum_s += wall_dt;

            constexpr int kMaxSubstepsPerFrame = 20;
            int substeps = 0;
            dropped_accum = false;

            while (accum_s >= dt && substeps < kMaxSubstepsPerFrame && !sim.isConcluded()) {
                sim.step(dt);
                refresh_obs();
                push_sample(simTime, last_obs);
                accum_s -= dt;
                ++substeps;
                advanced_this_frame = true;
            }

            last_substeps = substeps;

            if (substeps == kMaxSubstepsPerFrame) {
                accum_s = 0.0;
                dropped_accum = true;
            }
        } else {
            last_substeps = 0;
            dropped_accum = false;
        }

        if (!advanced_this_frame) {
            refresh_obs();
        }

        draft_vel_mps = v3((float)last_obs.draft_vel_mps_x,
                           (float)last_obs.draft_vel_mps_y,
                           (float)last_obs.draft_vel_mps_z);

        // --- ImGui frame ---
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

#ifndef VFEP_NO_IMGUI_DOCKING
        ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());
#endif

        // Enhanced Dashboard Overlay (Terminal/CMD Style - Professional)
        if (ui.show_hud) {
            ImGuiWindowFlags dashboard_flags =
                ImGuiWindowFlags_NoDecoration |
                ImGuiWindowFlags_AlwaysAutoResize |
                ImGuiWindowFlags_NoSavedSettings |
                ImGuiWindowFlags_NoFocusOnAppearing |
                ImGuiWindowFlags_NoNav;

            // Position at top-right
            ImVec2 viewport_size = ImGui::GetMainViewport()->Size;
            ImGui::SetNextWindowPos(ImVec2(viewport_size.x - 420, 12), ImGuiCond_Always, ImVec2(1.0f, 0.0f));
            ImGui::SetNextWindowBgAlpha(0.85f);

            if (ImGui::Begin("##Dashboard", &ui.show_hud, dashboard_flags)) {
                ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
                
                // Terminal green on dark background
                const ImVec4 header_col = ImVec4(0.0f, 1.0f, 0.0f, 1.0f);      // Bright green
                const ImVec4 metric_col = ImVec4(0.2f, 1.0f, 0.2f, 1.0f);      // Lighter green
                const ImVec4 status_ok = ImVec4(0.0f, 1.0f, 0.0f, 1.0f);       // Green
                const ImVec4 status_warn = ImVec4(1.0f, 1.0f, 0.0f, 1.0f);     // Yellow
                const ImVec4 status_fail = ImVec4(1.0f, 0.2f, 0.2f, 1.0f);     // Red
                
                // Title bar with branding
                ImGui::TextColored(header_col, "[ VFEP AUTONOMOUS SUPPRESSION ]");
                ImGui::SameLine(200);
                ImGui::TextColored(ImVec4(0.7f, 0.85f, 0.9f, 1.0f), "v1.0");
                ImGui::Separator();
                
                // Time + Status
                const bool concluded = sim.isConcluded();
                const char* state_txt = concluded ? "CONCLUDED" : (running ? "RUNNING" : "PAUSED");
                ImVec4 state_color = concluded ? status_fail : (running ? status_ok : ImVec4(0.5f, 0.8f, 1.0f, 1.0f));
                
                ImGui::Text("TIME: %.2f s", simTime);
                ImGui::SameLine(180);
                ImGui::TextColored(state_color, "[%s]", state_txt);
                ImGui::Spacing();
                
                // Fire Dynamics Section
                ImGui::TextColored(header_col, "=== FIRE DYNAMICS ===");
                
                double temp_C = (double)last_obs.T_K - 273.15;
                float temp_ratio = (float)std::clamp(temp_C / 600.0, 0.0, 1.0);
                ImGui::Text("Temp:  %.1f C", temp_C);
                ImGui::SameLine(180);
                ImGui::ProgressBar(temp_ratio, ImVec2(200, 12), "");
                
                double hrr_kW = 1e-3 * (double)last_obs.effective_HRR_W;
                float hrr_ratio = (float)std::clamp(hrr_kW / 500.0, 0.0, 1.0);
                ImGui::Text("HRR:   %.1f kW", hrr_kW);
                ImGui::SameLine(180);
                ImGui::ProgressBar(hrr_ratio, ImVec2(200, 12), "");
                ImGui::Spacing();
                
                // Suppression Performance Section
                ImGui::TextColored(header_col, "=== SUPPRESSION PERF ===");
                
                double kd = std::clamp((double)last_obs.knockdown_0_1, 0.0, 1.0);
                ImVec4 kd_color = (kd > 0.8) ? status_ok : (kd > 0.4) ? status_warn : status_fail;
                ImGui::Text("Knockdown: %.2f", kd);
                ImGui::SameLine(180);
                ImGui::ProgressBar((float)kd, ImVec2(200, 12), "");
                
                double hit_eff = std::clamp((double)last_obs.hit_efficiency_0_1, 0.0, 1.0);
                ImVec4 hit_color = (hit_eff > 0.8) ? status_ok : (hit_eff > 0.4) ? status_warn : status_fail;
                ImGui::Text("Hit Eff:    %.2f", hit_eff);
                ImGui::SameLine(180);
                ImGui::ProgressBar((float)hit_eff, ImVec2(200, 12), "");
                
                ImGui::Text("Agent flow: %.4f kg/s", (double)last_obs.agent_mdot_kgps);
                ImGui::Spacing();
                
                // Sector Coverage
                ImGui::TextColored(header_col, "=== SECTOR STATUS ===");
                int active_sectors = 0;
                for (int i = 0; i < 4; ++i) {
                    if (last_obs.sector_delivered_mdot_kgps[i] > 1e-6) active_sectors++;
                }
                ImGui::Text("Active: %d/4 sectors", active_sectors);
                
                // Draw mini sector boxes
                ImDrawList* draw_list = ImGui::GetWindowDrawList();
                ImVec2 sector_base = ImGui::GetCursorScreenPos();
                const float sector_size = 35.0f;
                const float sector_spacing = 50.0f;
                
                for (int i = 0; i < 4; ++i) {
                    double sector_kd = std::clamp(last_obs.sector_knockdown_0_1[i], 0.0, 1.0);
                    ImVec4 sector_color;
                    if (sector_kd > 0.7) sector_color = status_ok;
                    else if (sector_kd > 0.3) sector_color = status_warn;
                    else sector_color = status_fail;
                    
                    ImVec2 box_min = ImVec2(sector_base.x + i * sector_spacing, sector_base.y);
                    ImVec2 box_max = ImVec2(box_min.x + sector_size, box_min.y + sector_size);
                    ImU32 col = ImGui::GetColorU32(sector_color);
                    draw_list->AddRect(box_min, box_max, col, 2.0f, 0, 2.0f);
                    char sector_label[16];
                    snprintf(sector_label, sizeof(sector_label), "S%d", i);
                    draw_list->AddText(ImGui::GetIO().Fonts->Fonts[0], ImGui::GetFontSize(),
                                      ImVec2(box_min.x + 8, box_min.y + 12), col, sector_label);
                }
                ImGui::Dummy(ImVec2(sector_spacing * 4, sector_size + 10));
                ImGui::Spacing();
                
                // Regime indicator (colorized)
                ImGui::TextColored(header_col, "=== REGIME ===");
                const char* regime_text = suppression_regime_text(last_obs.suppression_regime);
                ImVec4 regime_color = status_ok;
                if (last_obs.suppression_regime == 1) regime_color = status_warn;
                else if (last_obs.suppression_regime == 0) regime_color = ImVec4(0.7f, 0.7f, 0.7f, 1.0f);
                
                ImGui::TextColored(regime_color, "[ %s ]", regime_text);
                
                if (dropped_accum) {
                    ImGui::Separator();
                    ImGui::TextColored(status_fail, ">> REALTIME DROPPED");
                }
                
                ImGui::PopStyleVar();
                ImGui::End();
            }
        }

        // Enhanced Control Console with Tabbed Interface
        if (ui.show_controls) {
            ImGui::SetNextWindowSize(ImVec2(600, 800), ImGuiCond_FirstUseEver);
            ImGui::Begin(">> CONTROL CONSOLE", &ui.show_controls);
            
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.2f, 1.0f, 0.2f, 1.0f));  // Terminal green
            ImGui::PushStyleColor(ImGuiCol_FrameBg, ImVec4(0.05f, 0.05f, 0.05f, 1.0f));
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.1f, 0.3f, 0.1f, 1.0f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.2f, 0.6f, 0.2f, 1.0f));
            
            const ImVec4 cmd_header = ImVec4(0.0f, 1.0f, 0.0f, 1.0f);
            const ImVec4 cmd_subsec = ImVec4(0.0f, 0.8f, 1.0f, 1.0f);
            
            // === TABBED INTERFACE ===
            if (ImGui::BeginTabBar("ControlTabs", ImGuiTabBarFlags_None)) {
                
                // ===== TAB 1: EXECUTION =====
                if (ImGui::BeginTabItem("  EXEC  ")) {
                    ImGui::TextColored(cmd_header, "[EXEC] Transport Controls");
                    ImGui::Separator();
                    
                    if (ImGui::Button(running ? "  PAUSE  " : "   RUN   ", ImVec2(100, 0))) running = !running;
                    ImGui::SameLine();
                    if (ImGui::Button("  STEP  ", ImVec2(100, 0))) {
                        if (!sim.isConcluded()) {
                            sim.step(dt);
                            refresh_obs();
                            push_sample(simTime, last_obs);
                            accum_s = 0.0;
                            last_substeps = 1;
                            dropped_accum = false;
                        }
                    }
                    ImGui::SameLine();
                    if (ImGui::Button(" RESET ", ImVec2(100, 0))) {
                        sim.resetToDataCenterRackScenario();
                        running = false;
                        refresh_obs();
                        accum_s = 0.0;
                        t_hist.clear(); T_hist.clear(); HRR_hist.clear(); O2_hist.clear();
                        EffExp_hist.clear(); KD_hist.clear(); KDTarget_hist.clear();
                        push_sample(simTime, last_obs);
                        last_substeps = 0;
                        dropped_accum = false;
                    }
                    
                    ImGui::Spacing();
                    float dt_slider = (float)dt;
                    ImGui::SliderFloat("Speed (1.0x)", &dt_slider, 0.005f, 0.200f, "%.3f s");
                    dt = (double)dt_slider;
                    ImGui::Spacing();
                    
                    // Scenario Selection
                    ImGui::TextColored(cmd_header, "[SCENARIO] Load Configuration");
                    ImGui::Separator();
                    
                    static int scenario_idx = 0;
                    static int agent_idx = 0;
                    
                    const char* scenario_names[] = {"Direct vs Glance", "Occlusion Wall", "Shielding Stack", "Mixed"};
                    ImGui::Combo(">> Scenario", &scenario_idx, scenario_names, IM_ARRAYSIZE(scenario_names));
                    
                    const char* agent_names[] = {"Clean Agent", "Dry Chemical", "CO2-like"};
                    ImGui::Combo(">> Agent", &agent_idx, agent_names, IM_ARRAYSIZE(agent_names));
                    
                    if (ImGui::Button("[ LOAD ]", ImVec2(-1, 0))) {
                        sim.resetToScenario((vfep::DemoScenario)scenario_idx, (vfep::AgentType)agent_idx);
                        running = false;
                        refresh_obs();
                        accum_s = 0.0;
                        t_hist.clear(); T_hist.clear(); HRR_hist.clear(); O2_hist.clear();
                        EffExp_hist.clear(); KD_hist.clear(); KDTarget_hist.clear();
                        push_sample(simTime, last_obs);
                        last_substeps = 0;
                        dropped_accum = false;
                    }
                    ImGui::Spacing();
                    
                    // Commands
                    ImGui::TextColored(cmd_header, "[COMMAND] System Actions");
                    ImGui::Separator();
                    
                    if (ImGui::Button("[ IGNITE ]", ImVec2(-1, 25))) {
                        if (!sim.isConcluded()) sim.commandIgniteOrIncreasePyrolysis();
                    }
                    if (ImGui::Button("[ START SUPPRESSION ]", ImVec2(-1, 25))) {
                        if (!sim.isConcluded()) sim.commandStartSuppression();
                    }
                    ImGui::Spacing();
                    
                    // Status Summary
                    ImGui::TextColored(cmd_header, "[STATUS] Current State");
                    ImGui::Separator();
                    
                    ImGui::Text("Time:        %.2f s", simTime);
                    ImGui::Text("Regime:      %s", suppression_regime_text(last_obs.suppression_regime));
                    ImGui::Text("HRR eff:     %.1f kW", 1e-3 * (double)last_obs.effective_HRR_W);
                    ImGui::Text("Knockdown:   %.1f %%", 100.0 * (double)last_obs.knockdown_0_1);
                    ImGui::Text("Hit eff:     %.1f %%", 100.0 * (double)last_obs.hit_efficiency_0_1);
                    ImGui::Text("Agent flow:  %.4f kg/s", (double)last_obs.agent_mdot_kgps);
                    
                    if (last_substeps > 0) {
                        ImGui::TextColored(ImVec4(0.0f, 1.0f, 1.0f, 1.0f), "Substeps:    %d", last_substeps);
                    }
                    
                    ImGui::EndTabItem();
                }
                
                // ===== TAB 2: NOZZLE =====
                if (ImGui::BeginTabItem("  NOZZLE  ")) {
                    ImGui::TextColored(cmd_header, "[NOZZLE] Pose Control");
                    ImGui::Separator();
                    
                    // Manual pose controls only available when override is OFF (not using rail)
                    if (viz_override_nozzle_pose) {
                        ImGui::TextDisabled("Position locked to rail (s=%.2f)", viz_nozzle_s_0_1);
                        ImGui::TextDisabled("Use sliders below to control nozzle");
                    } else {
                        ImGui::DragFloat3("Position (m)##noz", &nozzle_pos.x, 0.05f);
                        ImGui::DragFloat3("Direction##noz", &nozzle_dir.x, 0.02f);
                        
                        if (ImGui::Button("[ APPLY NOZZLE ]", ImVec2(-1, 0))) {
                            sim.setNozzlePose({(double)nozzle_pos.x, (double)nozzle_pos.y, (double)nozzle_pos.z},
                                              {(double)nozzle_dir.x, (double)nozzle_dir.y, (double)nozzle_dir.z});
                            refresh_obs();
                        }
                    }
                    ImGui::Spacing();
                    
                    ImGui::TextColored(cmd_header, "[DEBUG] Rail Parameters");
                    ImGui::Separator();
                    ImGui::DragFloat("Rail drop (m)", &rail_ceiling_drop_m, 0.01f, 0.0f, 2.0f, "%.2f");
                    ImGui::DragFloat("Rail margin (m)", &rail_margin_m, 0.01f, 0.0f, 2.0f, "%.2f");
                    ImGui::DragFloat("Nozzle drop (m)", &nozzle_drop_from_rail_m, 0.01f, 0.0f, 2.0f, "%.2f");
                    
                    ImGui::Checkbox("Override nozzle pose", &viz_override_nozzle_pose);
                    if (viz_override_nozzle_pose) {
                        ImGui::SliderFloat("Nozzle s (0-1)", &viz_nozzle_s_0_1, 0.0f, 1.0f, "%.3f");
                        ImGui::SliderFloat("Nozzle pan (deg)", &viz_nozzle_pan_deg, -180.0f, 180.0f, "%.1f");
                        ImGui::SliderFloat("Nozzle tilt (deg)", &viz_nozzle_tilt_deg, -90.0f, 90.0f, "%.1f");
                    }
                    
                    ImGui::Text("Rail valid: %s", ceiling_rail.isValid() ? "YES" : "NO");
                    ImGui::Text("Nozzle valid: %s", rail_nozzle.isValid() ? "YES" : "NO");
                    
                    ImGui::EndTabItem();
                }
                
                // ===== TAB 3: VISUALIZATION =====
                if (ImGui::BeginTabItem("  VIZ  ")) {
                    ImGui::TextColored(cmd_header, "[VISUALIZATION] Draw Layers");
                    ImGui::Separator();
                    
                    ImGui::Checkbox("Warehouse", &ui.draw_warehouse);
                    ImGui::Checkbox("Rack", &ui.draw_rack);
                    ImGui::Checkbox("Fire volume", &ui.draw_fire);
                    ImGui::Checkbox("Fire sectors", &ui.draw_fire_sectors);
                    ImGui::Checkbox("Ceiling Rail", &ui.draw_ceiling_rail);
                    ImGui::Checkbox("Nozzle marker", &ui.draw_nozzle);
                    ImGui::Checkbox("Spray cone", &ui.draw_spray);
                    ImGui::Checkbox("Hit marker", &ui.draw_hit_marker);
                    ImGui::Checkbox("Draft arrow", &ui.draw_draft);
                    
                    ImGui::EndTabItem();
                }
                
                // ===== TAB 4: PLOTS =====
                if (ImGui::BeginTabItem("  PLOTS  ")) {
                    const int N = (int)t_hist.size();
                    const int start = (N > kPlotWindowN) ? (N - kPlotWindowN) : 0;
                    const int count = N - start;

                    if (count > 1) {
                        const double t0 = t_hist[start];
                        const double t1 = t_hist[start + count - 1];
                        ImGui::Text("Samples: %d   Window: [%0.2f, %0.2f] s", N, t0, t1);
                        ImGui::Separator();

                        plot_line_with_xlimits("Temperature (K)", "T_K",
                                               t_hist.data() + start, T_hist.data() + start, count, t0, t1);

                        plot_line_with_xlimits("HRR (W)", "HRR_W",
                                               t_hist.data() + start, HRR_hist.data() + start, count, t0, t1);

                        plot_line_with_xlimits("Effective Exposure (kg)", "EffExp_kg",
                                               t_hist.data() + start, EffExp_hist.data() + start, count, t0, t1);

                        if (ImPlot::BeginPlot("Knockdown (0-1)")) {
                            #if defined(ImAxis_X1)
                            ImPlot::SetupAxisLimits(ImAxis_X1, t0, t1, ImGuiCond_Always);
                            #elif defined(ImPlotAxis_X1)
                            ImPlot::SetupAxisLimits(ImPlotAxis_X1, t0, t1, ImGuiCond_Always);
                            #else
                            #endif

                            ImPlot::PlotLine("KD", t_hist.data() + start, KD_hist.data() + start, count);
                            ImPlot::PlotLine("KD_target", t_hist.data() + start, KDTarget_hist.data() + start, count);

                            ImPlot::EndPlot();
                        }

                        plot_line_with_xlimits("O2 (vol %)", "O2",
                                               t_hist.data() + start, O2_hist.data() + start, count, t0, t1);
                    } else {
                        ImGui::Text("Samples: %d", N);
                        ImGui::TextUnformatted("No data yet (press Run or Step).");
                    }
                    
                    ImGui::EndTabItem();
                }
                
                ImGui::EndTabBar();
            }
            
            ImGui::Spacing();
            ImGui::Separator();
            ImGui::TextColored(ImVec4(0.5f, 0.7f, 0.5f, 1.0f), ">> VFEP Simulation v1.0");
            ImGui::TextColored(ImVec4(0.4f, 0.6f, 0.4f, 1.0f), ">> Microgravity Fire Suppression");
            
            ImGui::PopStyleColor(4);
            ImGui::End();
        }

        // --- Step 2: model-backed rail-mounted nozzle (viz-driven, kinematics only) ---
        // Recompute the rail every frame because the nozzle pose depends on it.
        // (We will draw the rail later; this block is about producing pose.)
        ceiling_rail_cfg.drop_from_ceiling_m = (double)rail_ceiling_drop_m;
        ceiling_rail_cfg.margin_from_rack_m  = (double)rail_margin_m;
        ceiling_rail.setConfig(ceiling_rail_cfg);

        ceiling_rail_in.ceiling_y_m      = 0.0; // keep default behavior (ceiling_y = 2*warehouse_half.y)
        ceiling_rail_in.warehouse_half_m = to_v3d(warehouse_half);
        ceiling_rail_in.rack_center_m    = to_v3d(rack_center);
        ceiling_rail_in.rack_half_m      = to_v3d(rack_half);
        ceiling_rail.recompute(ceiling_rail_in);

        // If UI override is enabled, compute nozzle pose from rail + s + yaw/pitch + drop.
        if (viz_override_nozzle_pose && ceiling_rail.isValid()) {
            rail_nozzle_cfg.nozzle_drop_from_rail_m = (double)nozzle_drop_from_rail_m;
            rail_nozzle.setConfig(rail_nozzle_cfg);

            rail_nozzle_in.ceiling_rail = &ceiling_rail;
            rail_nozzle_in.s_0_1        = (double)viz_nozzle_s_0_1;
            rail_nozzle_in.yaw_deg      = (double)viz_nozzle_pan_deg;
            rail_nozzle_in.pitch_deg    = (double)viz_nozzle_tilt_deg;

            rail_nozzle.recompute(rail_nozzle_in);

            if (rail_nozzle.isValid()) {
                nozzle_pos = to_v3f(rail_nozzle.pose().nozzle_pos_room_m);
                nozzle_dir = to_v3f(rail_nozzle.pose().spray_dir_unit_room);
            }
        } else if (!viz_override_nozzle_pose && last_obs.agent_mdot_kgps > 1e-6) {
            // When not overriding and suppression is active, aim toward hotspot
            const auto np = sim.getNozzlePos_m();
            nozzle_pos = v3((float)np.x, (float)np.y, (float)np.z);
            
            // Auto-aim: compute direction from nozzle to fire hotspot
            Vec3f to_fire = sub(fire_center, nozzle_pos);
            const float dist = len(to_fire);
            if (dist > 1e-3f) {
                nozzle_dir = mul(to_fire, 1.0f / dist);  // normalize
                // Update simulation with auto-aim direction
                sim.setNozzlePose({(double)nozzle_pos.x, (double)nozzle_pos.y, (double)nozzle_pos.z},
                                  {(double)nozzle_dir.x, (double)nozzle_dir.y, (double)nozzle_dir.z});
            } else {
                // Fallback to observed spray direction if fire too close
                nozzle_dir = v3((float)last_obs.spray_dir_unit_x,
                               (float)last_obs.spray_dir_unit_y,
                               (float)last_obs.spray_dir_unit_z);
            }
        }

        ImGui::Render();

        int fb_w = 0, fb_h = 0;
        glfwGetFramebufferSize(window, &fb_w, &fb_h);

        if (fb_w > 0 && fb_h > 0) {
            glViewport(0, 0, fb_w, fb_h);

            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_LESS);
            glDisable(GL_CULL_FACE);

            glClearColor(0.06f, 0.06f, 0.07f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            const float aspect = (fb_h > 0) ? (float)fb_w / (float)fb_h : 1.0f;
            set_perspective(55.0f, aspect, 0.05f, 100.0f);

            const float yaw   = cam_yaw_deg   * 3.1415926535f / 180.0f;
            const float pitch = cam_pitch_deg * 3.1415926535f / 180.0f;

            Vec3f eye = v3(
                cam_target.x + cam_dist * std::cos(pitch) * std::sin(yaw),
                cam_target.y + cam_dist * std::sin(pitch),
                cam_target.z + cam_dist * std::cos(pitch) * std::cos(yaw)
            );
            look_at(eye, cam_target, v3(0.0f, 1.0f, 0.0f));

            if (ui.draw_warehouse) {
                glColor3f(0.25f, 0.25f, 0.28f);
                draw_wire_box(v3(0.0f, warehouse_half.y, 0.0f), warehouse_half);
            }

            if (ui.draw_rack) {
                const float rackTempC = (float)(last_obs.T_K - 273.15);
                float rr, rg, rb;
                temp_to_color(rackTempC, rr, rg, rb);
                glColor3f(rr, rg, rb);
                draw_solid_box(rack_center, rack_half);

                glColor3f(0.05f, 0.05f, 0.05f);
                draw_wire_box(rack_center, rack_half);
            }

            if (ui.draw_ceiling_rail) {
                // Push UI params into model config
                ceiling_rail_cfg.drop_from_ceiling_m = (double)rail_ceiling_drop_m;
                ceiling_rail_cfg.margin_from_rack_m  = (double)rail_margin_m;
                ceiling_rail.setConfig(ceiling_rail_cfg);

                // Push scene geometry into model inputs.
                // NOTE: leaving ceiling_y_m = 0 uses default: ceiling_y = 2*warehouse_half.y
                //       which matches the legacy draw_square_rail() behavior.
                ceiling_rail_in.ceiling_y_m      = 0.0;
                ceiling_rail_in.warehouse_half_m = to_v3d(warehouse_half);
                ceiling_rail_in.rack_center_m    = to_v3d(rack_center);
                ceiling_rail_in.rack_half_m      = to_v3d(rack_half);

                ceiling_rail.recompute(ceiling_rail_in);

                if (ceiling_rail.isValid()) {
                    glColor3f(0.85f, 0.85f, 0.15f);

                    const auto& g = ceiling_rail.geometry();
                    const Vec3f p00 = to_v3f(g.corners_room_m[0]);
                    const Vec3f p10 = to_v3f(g.corners_room_m[1]);
                    const Vec3f p11 = to_v3f(g.corners_room_m[2]);
                    const Vec3f p01 = to_v3f(g.corners_room_m[3]);

                    draw_line(p00, p10);
                    draw_line(p10, p11);
                    draw_line(p11, p01);
                    draw_line(p01, p00);
                }
            }

            if (ui.draw_fire) {
                const double hrr_vis_W = (std::isfinite(last_obs.effective_HRR_W) && last_obs.effective_HRR_W > 0.0)
                    ? last_obs.effective_HRR_W
                    : last_obs.HRR_W;

                const float fire_s = fire_scale_from_HRR_W(hrr_vis_W);
                Vec3f fire_half = mul(v3(0.35f, 0.45f, 0.35f), fire_s);

                if (hrr_vis_W > 1.0) {
                    if (ui.draw_fire_sectors) {
                        Vec3f sub_half = v3(fire_half.x * 0.48f, fire_half.y, fire_half.z * 0.48f);

                        const float sx[4] = {-1.0f, +1.0f, -1.0f, +1.0f};
                        const float sz[4] = {-1.0f, -1.0f, +1.0f, +1.0f};

                        for (int i = 0; i < 4; ++i) {
                            const float kd = clampf((float)last_obs.sector_knockdown_0_1[i], 0.0f, 1.0f);
                            const float intensity = clampf(0.20f + 0.80f * (1.0f - kd), 0.0f, 1.0f);

                            glColor3f(0.85f * intensity, 0.25f * intensity, 0.05f * intensity);

                            Vec3f c = fire_center;
                            c.x += sx[i] * sub_half.x;
                            c.z += sz[i] * sub_half.z;
                            draw_solid_box(c, sub_half);
                        }
                    }

                    glColor3f(0.15f, 0.05f, 0.02f);
                    draw_wire_box(fire_center, fire_half);
                }
            }

            if (ui.draw_draft) {
                const float mag = len(draft_vel_mps);
                const float L = clampf(draft_arrow_scale * mag, 0.2f, 4.0f);
                glColor3f(0.10f, 0.75f, 0.75f);
                draw_arrow(add(rack_center, v3(0.0f, rack_half.y + 0.3f, 0.0f)), draft_vel_mps, L);
            }

            if (ui.draw_nozzle) {
                glColor3f(0.55f, 0.55f, 0.60f);
                draw_wire_box(nozzle_pos, v3(0.06f, 0.06f, 0.06f));
            }

            const float eff_draw = clampf((float)last_obs.hit_efficiency_0_1, 0.0f, 1.0f);
            const float cone_len = clampf(spray_L0 + spray_L1 * eff_draw, 0.0f, spray_max_len);
            const float cone_rad = clampf(spray_R0 + spray_R1 * eff_draw, 0.0f, 3.0f);

            const Vec3f eff_dir = norm(v3((float)last_obs.spray_dir_unit_x,
                                          (float)last_obs.spray_dir_unit_y,
                                          (float)last_obs.spray_dir_unit_z));

            const Vec3f nozzle_dir_n = (len(nozzle_dir) > 1e-6f) ? norm(nozzle_dir) : eff_dir;

            if (ui.draw_spray && last_obs.agent_mdot_kgps > 1e-6) {
                glColor3f(0.55f, 0.25f, 0.70f);
                draw_cone_world(nozzle_pos, eff_dir, cone_len, cone_rad, 18);

                glColor3f(0.35f, 0.18f, 0.45f);
                draw_line(nozzle_pos, add(nozzle_pos, mul(eff_dir, cone_len)));

                glColor3f(0.22f, 0.22f, 0.25f);
                draw_line(nozzle_pos, add(nozzle_pos, mul(nozzle_dir_n, cone_len)));
            }

            if (ui.draw_hit_marker) {
                float t_hit = 0.0f;
                if (len(eff_dir) > 1e-6f && ray_aabb_intersect(nozzle_pos, eff_dir, rack_center, rack_half, t_hit)) {
                    Vec3f hit = add(nozzle_pos, mul(eff_dir, t_hit));

                    const float marker_half = clampf(hit_marker_base + hit_marker_gain * eff_draw, 0.01f, 0.5f);
                    Vec3f mh = v3(marker_half, marker_half, marker_half);

                    const float hit_quality = eff_draw;
                    const float r = 0.80f * (1.0f - hit_quality) + 0.25f * hit_quality;
                    const float g = 0.20f * (1.0f - hit_quality) + 0.90f * hit_quality;
                    const float b = 0.75f * (1.0f - hit_quality) + 0.30f * hit_quality;

                    glColor3f(r, g, b);
                    draw_solid_box(hit, mh);

                    glColor3f(0.08f, 0.03f, 0.10f);
                    draw_wire_box(hit, mh);
                }
            }

            glDisable(GL_DEPTH_TEST);
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        }

        glfwSwapBuffers(window);
    }

    // Cleanup
    if (implot_ctx) ImPlot::DestroyContext();
    if (imgui_gl3) ImGui_ImplOpenGL3_Shutdown();
    if (imgui_glfw) ImGui_ImplGlfw_Shutdown();
    if (imgui_ctx) ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
