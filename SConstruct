import os
env = Environment(ENV = os.environ)

env.Append(LIBS=["png"])
env.Append(CXXFLAGS=["-std=c++11","-g","-Wall","-O3","-I/opt/homebrew/Cellar/libpng/1.6.37/include"])
env.Append(LINKFLAGS=["-L/opt/homebrew/Cellar/libpng/1.6.37/lib"])


env.Program("driver",["main.cpp","parse.cpp","dump_png.cpp","driver_state.cpp","shaders.cpp"])
