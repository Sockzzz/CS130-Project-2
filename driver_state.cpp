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
    //fixed initialization
    state.image_depth = new float[width*height];
    state.image_color = new pixel[width*height];
    for(int i = 0; i < width * height; ++i){
        state.image_color[i] = make_pixel(0,0,0);
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
    for(int i = 0; i < startingPoints.size(); ++i){

        data_geometry toAdd;
        toAdd.data = startingPoints.at(i).data;

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


data_geometry giveIntersection(data_geometry a, data_geometry b){

    //all the points needed for the math
    float v0_x = a.gl_Position[0];
    float v0_y = a.gl_Position[1];
    float v1_x = b.gl_Position[0];
    float v1_y = b.gl_Position[1];

    //equation of a line stuff
    float rise = v1_y - v0_y;
    float run = v1_x - v0_x;
    float slope = rise/run;

    //y intercept
    float intercept = v1_y - (slope*v1_x);

    //4 sub cases, depending on what way





}


void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
/*
    //data_geometry has
        //glposition vec4
        //data pointer

    //case 1: everything is inside
    //case 2: one vertex is outside
    //case 3: 2 vertices are outside
    //case 4: all vertices are outside
        //4.1: entire triangle is outside
        //4.2: triangle is still in viewing area but vertices are outside

    //check to see how many are inside
    //range should be -1 -> 1 for both vertices
    int needs_correction;

    //easier handling for the geometry
    vector<data_geometry> vertices;
    vertices.push_back(v0);
    vertices.push_back(v1);
    vertices.push_back(v2);


    //stores which DG are outside
    bool OOB[3];

    //stores whether the x or y or both are OOB, will have 3 bool[2]s
    vector<bool[2]> control;

    //used for laziness
    int outside = 0;

    //determining how much is out of bounds
    for(int i = 0; i < 3; i++){
        bool tosend[2];
        control.push_back(tosend);

        for(int j = 0; j < 2; ++ j){
            float temp = vertices.at(i).gl_Position[j];

            if(temp > 1 || temp < -1) {
                OOB[i] = true;
                control.at(i)[j] = true;
            }
            else
                OOB[i] = false;
                control.at(i)[j] = false;
        }
    }

    //counts how many DG will need to be fixed
    for(int b = 0; b < 3; ++b){
        if(OOB[b])
            ++outside;
    }

    if(outside == 0){
        face = 6;
    }
    else if(outside == 1){

        //determine which DG is outside
        //determine if x or/and y is outside
        //then reroute

        int whichDG;
        for(int i = 0 ; i < 3; ++i){
            if(OOB[i])
                whichDG = i;
        }

        //use midpoint and get the new coordinates
        //send new coordinates + (v[i] + 1) to clip again


    }


    ////////////////////////////////////////////////////////////
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,v0,v1,v2,face+1);
    */
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{

    int h = state.image_height/2;
    int w = state.image_width/2;



    //have to scale the coordinates correctly to get the in NCD
    int x0 = (v0.gl_Position[0]/v0.gl_Position[3]* w) + w;
    int y0 = (v0.gl_Position[1]/v0.gl_Position[3]* h) + h;

    int x1 = (v1.gl_Position[0]/v1.gl_Position[3] * w) + w;
    int y1 = (v1.gl_Position[1]/v1.gl_Position[3] * h) + h;

    int x2 = (v2.gl_Position[0]/v2.gl_Position[3] * w) + w;
    int y2 = (v2.gl_Position[1]/v2.gl_Position[3] * h) + h;



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


            //changed the [2] from magnitude because we want to get the Z axis since we need that for the interpolation
            double area_total = cross(vertex_2-vertex_1, vertex_2-vertex_3)[2] * .5;
            double area_pbc = cross(vertex_2 - point,vertex_2-vertex_3)[2] * .5;
            double area_apc = cross(point-vertex_1,point-vertex_3)[2] * .5;
            double area_abp = cross(vertex_2-vertex_1,vertex_2-point)[2] * .5;

            double alpha = area_pbc / area_total;
            double beta = area_apc / area_total;
            double gamma = area_abp / area_total;

            //the z buffer point scalers
            float scale_a = v0.gl_Position[2]/v0.gl_Position[3];
            float scale_b = v1.gl_Position[2]/v1.gl_Position[3];
            float scale_c = v2.gl_Position[2]/v2.gl_Position[3];

            //this is the points true z cordinate that we can then compare to for zbuff
            point[2] = alpha*scale_a + beta*scale_b + gamma*scale_c;

            //bounds check
            if(alpha >= 0 && beta >= 0 && gamma >= 0){

                //gonna get sent to the fragment shader for
                data_fragment topass;
                topass.data = new float[MAX_FLOATS_PER_VERTEX];


                //determine and set frag type
                for(int a = 0; a < MAX_FLOATS_PER_VERTEX; ++a){

                    interp_type chosenOne = state.interp_rules[a];

                    //flat shader so we just give it the first color and it copies
                    //needs [a] or doesnt work for test 8, idk why
                    if(chosenOne == interp_type::flat){
                        topass.data[a] = v0.data[a];
                    }

                    //incomplete
                    else if(chosenOne == interp_type::smooth){
                        double scalar = alpha/v0.gl_Position[3] + beta/v1.gl_Position[3] + gamma/v2.gl_Position[3];

                        double alpha_s = alpha/(scalar*v0.gl_Position[3]);
                        double beta_s = beta/(scalar*v1.gl_Position[3]);
                        double gamma_s = gamma/(scalar*v2.gl_Position[3]);

                        topass.data[a] = alpha_s*v0.data[a] + beta_s*v1.data[a] + gamma_s*v2.data[a];

                    }
                    //incomplete
                    else if(chosenOne == interp_type::invalid){
                        //need to implement invalid later not there yet
                    }
                    //know this, need to give it the barry affected color coordinates, once for each color channel
                    else if(chosenOne == interp_type::noperspective){
                        topass.data[a] = alpha*v0.data[a] + beta*v1.data[a] + gamma*v2.data[a];
                    }
                    //just in case something unintended happens
                    else{
                        exit(309);
                    }

                }


                // a = (b/g + c/g + d/g)*g == b + c + d

                //create this to throw it into frag shader
                data_output out;
                state.fragment_shader(topass,out,state.uniform_data);

                //compare z cordinate to buffer, if smaller than buffer, draw it and update, otherwise dont draw
                if(state.image_depth[x+(y*state.image_width)] > point[2]) {
                    state.image_depth[x+(y*state.image_width)] = point[2];
                    //cout<<"point[2]: "<<point[2]<<endl;
                    //cout<<"image depth after: "<<*state.image_depth<<endl;

                    state.image_color[x + (y * state.image_width)] = make_pixel(255 * out.output_color[0],
                                                                                255 * out.output_color[1],
                                                                                255 * out.output_color[2]);

                }


            }


        }
    }


}

