#include "driver_state.h"
#include <cstring>
#include <limits>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;

    int numEntries = width * height;
    state.image_color = new pixel[numEntries];
    state.image_depth = new float[numEntries];
    for (int i = 0; i < numEntries; i++) {
        pixel newPix = make_pixel(0,0,0);
        state.image_color[i] = newPix;
        state.image_depth[i] = std::numeric_limits<float>::max();
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    data_geometry * g_data = new data_geometry[state.num_vertices];
    for (int i = 0; i < state.num_vertices; i++) {
        g_data[i].gl_Position = vec4(0,0,0,0);
        g_data[i].data = new float(0.0);    // is this good enough lol
        
        data_vertex v_data; v_data.data = &state.vertex_data[state.floats_per_vertex * i];
        state.vertex_shader(v_data, g_data[i], state.uniform_data);
    }

    switch (type) {
        case render_type::triangle:
            for (int i = 0; i < (state.num_vertices/* * state.floats_per_vertex*/); i += 3) {
                rasterize_triangle(state, g_data[i], g_data[i+1], g_data[i+2]);
            }
            break;
        case render_type::indexed:
            break;
        case render_type::fan:
            break;
        case render_type::strip:
            break;
    }

    // delete [] g_data;
    for (int i = 0; i < state.num_vertices; i++) {
        delete g_data[i].data;
    }

    delete [] g_data;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,v0,v1,v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    // rasterizes a SINGLE triangle
    // std::cout<<"TODO: implement rasterization"<<std::endl;
    // Vertex Positions
    vec4 A(
                v0.gl_Position[0] / v0.gl_Position[3],  // X
                v0.gl_Position[1] / v0.gl_Position[3],  // Y
                v0.gl_Position[2] / v0.gl_Position[3],  // Z
                v0.gl_Position[3]                       // W
                );
    vec4 B(
                v1.gl_Position[0] / v1.gl_Position[3],
                v1.gl_Position[1] / v1.gl_Position[3],
                v1.gl_Position[2] / v1.gl_Position[3],
                v1.gl_Position[3]
                );
    vec4 C(
                v2.gl_Position[0] / v2.gl_Position[3],
                v2.gl_Position[1] / v2.gl_Position[3],
                v2.gl_Position[2] / v2.gl_Position[3],
                v2.gl_Position[3]
                );
    
    A[0] = (state.image_width - 1) * ((A[0] + 1) / 2.0);  // image_width -> indices, scale the x/y to be 0 to 1 (from -1 to 1)
    A[1] = (state.image_height - 1) * ((A[1] + 1) / 2.0);
    
    B[0] = (state.image_width - 1) * ((B[0] + 1) / 2.0);
    B[1] = (state.image_height - 1) * ((B[1] + 1) / 2.0);

    C[0] = (state.image_width - 1) * ((C[0] + 1) / 2.0);
    C[1] = (state.image_height - 1) * ((C[1] + 1) / 2.0);

    double alpha, beta, gamma;

    const double Tri_Area = 0.5 * (
                            (B[0] * C[1] - C[0] * B[1]) +   // B_x * C_y - C_x * B_y
                            (C[0] * A[1] - A[0] * C[1]) +   // C_x * A_y - A_x * C_y
                            (A[0] * B[1] - B[0] * A[1])     // A_x * B_y - B_x * A_y
                            );

    for (int i = 0; i < state.image_width; i++) {
        for (int j = 0; j < state.image_height; j++) {
            int pixel_pos = i + (j * state.image_width);

            alpha = 0.5 * (
                            (B[0] * C[1] - C[0] * B[1]) +   // B_x * C_y - C_x * B_y
                            (C[0] * j - i * C[1]) +         // C_x * P_y - P_x * C_y
                            (i * B[1] - B[0] * j)           // P_x * B_y - B_x * P_y
                            );
            beta = 0.5 * (
                            (i * C[1] - C[0] * j) +         // P_x * C_y - C_x * P_y
                            (C[0] * A[1] - A[0] * C[1]) +   // C_x * A_y - A_x * C_y
                            (A[0] * j - i * A[1])           // A_x * P_y - P_x * A_y
                            );
            gamma = 0.5 * (
                            (B[0] * j - i * B[1]) +         // B_x * P_y - P_x * B_y
                            (i * A[1] - A[0] * j) +         // P_x * A_y - A_x * P_y
                            (A[0] * B[1] - B[0] * A[1])     // A_x * B_y - B_x * A_y
                            );

            alpha /= Tri_Area;
            beta /= Tri_Area;
            gamma /= Tri_Area;

            if (alpha >= -0.001 && beta >= -0.001 && gamma >= -0.001) {
                data_fragment f_data;   // store fragment data (same as data_vertex, just let the function do the work)
                data_output out_data;   // vec4 of rgb (determined by fragment_shader)
                f_data.data = new float[MAX_FLOATS_PER_VERTEX];

                state.fragment_shader(f_data, out_data, state.uniform_data);
                float r = out_data.output_color[0], g = out_data.output_color[1], b = out_data.output_color[2];

                float temp_depth = alpha * A[2] + beta * B[2] + gamma * C[2];

                if (state.floats_per_vertex > 3) {
                    // Denomenator of eq. for perspective correct barycentric: a'/w_a + b'/w_b + c'/w_c
                    double bay_w = (alpha / A[3]) + (beta / B[3]) + (gamma / C[3]);  
                    
                    switch(state.interp_rules[3]) {  // for the color mixes
                        case interp_type::flat:
                            // test 07
                            r = v0.data[3];
                            g = v0.data[4];
                            b = v0.data[5];
                            break;
                        case interp_type::smooth:   // perspective correct interpolation. need k and w
                           
                            // Calculating correct alpha, beta, gamma: x'/w_x / (a'/w_a + b'/w_b + c'/w_c)
                            alpha = (alpha / A[3]) / bay_w;
                            beta = (beta / B[3]) / bay_w;
                            gamma = (gamma / C[3]) / bay_w;

                            // apply barycentric with perspective correct alpha, beta, gamma
                            r = v0.data[3] * alpha + v1.data[3] * beta + v2.data[3] * gamma;
                            g = v0.data[4] * alpha + v1.data[4] * beta + v2.data[4] * gamma;
                            b = v0.data[5] * alpha + v1.data[5] * beta + v2.data[5] * gamma;
                            break;
                        case interp_type::noperspective:    // nonperspective barycentric for interpolation
                            r = v0.data[3] * alpha + v1.data[3] * beta + v2.data[3] * gamma;
                            g = v0.data[4] * alpha + v1.data[4] * beta + v2.data[4] * gamma;
                            b = v0.data[5] * alpha + v1.data[5] * beta + v2.data[5] * gamma;
                            break;
                    }
                }
                if (temp_depth < state.image_depth[pixel_pos]) {
                    state.image_depth[pixel_pos] = temp_depth;
                    state.image_color[pixel_pos] = make_pixel(r * 255, g * 255, b * 255);
                    // state.image_color[i + (j * state.image_width)] = make_pixel(r, g, b);
                    // state.image_color[i + (j * state.image_width)] = make_pixel(out_data.output_color[0] * 255, 
                    //                                                             out_data.output_color[1] * 255, 
                    //                                                             out_data.output_color[2] * 255);
                }
            }
        }
    }
}

