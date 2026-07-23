cmake_minimum_required(VERSION 4.2)
project(sensor_fusion_playground C)

set(CMAKE_C_STANDARD 17)

add_executable(sensor_fusion_playground
        src/main.c
        src/path.c
        src/bspline.c
        src/matrix.c
)

target_include_directories(sensor_fusion_playground PRIVATE inc)
