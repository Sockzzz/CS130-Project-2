#include "driver_state.h"
#include <cstring>
#include <vector>

using namespace std;

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

    //image color is a pixel *
        //color data is stored
        //first image_width entries correspond to bottom row of the image
        //the next image_width entries correspond to the bottom row
        //array has image_width*image_height entries

        //pixel *
            //packs RGBA vales
            //use make/from pixel to pack or unpack data

    //image depth is a float *
        //the array stores the depth of a pixel and is used for z buffering
        //size and layout is the same as image_color

    //make the whole screen black
    state.image_color = new pixel[width*height];
    for(int i = 0; i < width * height; ++i){
        state.image_color[i] = make_pixel(0,0,0);
    }

    state.image_depth = new float[width*height];
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

    //set up structures to hold data_vertex and data_geometry
    vector<data_vertex> startingPoints;
    vector<data_geometry> vertices;
    int toiterate = state.num_vertices;

    //fill up the vertex objects
    for(int i = 0; i < toiterate; ++i){

        //slide selector to right spot
        int space = i * state.floats_per_vertex;

        //create object and add to stack
        data_vertex toAdd;
        toAdd.data = state.vertex_data + space;
        startingPoints.push_back(toAdd);

    }

    //then i call vertex shader to build data_geometry
    for(int i = 0; i < startingPoints.size(); ++ i){

        data_geometry toAdd;

        //call generic vertex shader, use 0
        state.vertex_shader(startingPoints.at(i),toAdd,state.uniform_data);
        vertices.push_back(toAdd);

    }


    //then i call rasterize, need to loop to account for multiple
    for(int i = 0; i <= vertices.size()-3; i+=3){

        rasterize_triangle(state,vertices.at(i),
                           vertices.at(i+1),vertices.at(i+2));

    }

}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.

//
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

    int h = state.image_height/2;
    int w = state.image_width/2;

    int x0 = (v0.gl_Position[0]* w) + w;
    int y0 = (v0.gl_Position[1]* h) + h;

    int x1 = (v1.gl_Position[0] * w) + w;
    int y1 = (v1.gl_Position[1] * h) + h;

    int x2 = (v2.gl_Position[0] * w) + w;
    int y2 = (v2.gl_Position[1] * h) + h;

    int min_x = fmin(x0, fmin(x1,x2));
    int max_x = fmax(x0, fmax(x1,x2));

    int min_y = fmin(y0, fmin(y1,y2));
    int max_y = fmax(y0, fmax(y1, y2));

    for(int x = min_x; x < max_x; ++x){
        for(int y = min_y; y < max_y; ++y){

            vec3 vertex_1 = vec3(x0,y0,0);
            vec3 vertex_2 = vec3(x1, y1, 0);
            vec3 vertex_3 = vec3(x2, y2, 0);
            vec3 point = vec3(x,y,0);

            double area_total = cross(vertex_2-vertex_1, vertex_2-vertex_3).magnitude() * .5;
            double area_pbc = cross(vertex_2 - point,vertex_2-vertex_3).magnitude() * .5;
            double area_apc = cross(point-vertex_1,point-vertex_3).magnitude() * .5;
            double area_abp = cross(vertex_2-vertex_1,vertex_2-point).magnitude() * .5;

            double alpha = area_pbc / area_total;
            double beta = area_apc / area_total;
            double gamma = area_abp / area_total;

            if(alpha + beta + gamma <= 1){

                state.image_color[x+(y*state.image_width)] = make_pixel(255,255,255);

            }


        }
    }


}

