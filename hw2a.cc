#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define PNG_NO_SETJMP
#include <sched.h>
#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <pthread.h>

pthread_mutex_t mutex;
int i_cur = 0, j_cur = 0;

void write_png(const char* filename, int iters, int width, int height, const int* buffer) {
    FILE* fp = fopen(filename, "wb");
    assert(fp);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    assert(png_ptr);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    assert(info_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_filter(png_ptr, 0, PNG_NO_FILTERS);
    png_write_info(png_ptr, info_ptr);
    png_set_compression_level(png_ptr, 1);
    size_t row_size = 3 * width * sizeof(png_byte);
    png_bytep row = (png_bytep)malloc(row_size);
    for (int y = 0; y < height; ++y) {
        memset(row, 0, row_size);
        for (int x = 0; x < width; ++x) {
            int p = buffer[(height - 1 - y) * width + x];
            png_bytep color = row + x * 3;
            if (p != iters) {
                if (p & 16) {
                    color[0] = 240;
                    color[1] = color[2] = p % 16 * 16;
                } else {
                    color[0] = p % 16 * 16;
                }
            }
        }
        png_write_row(png_ptr, row);
    }
    free(row);
    png_write_end(png_ptr, NULL);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
}

typedef struct{
    const char* filename;
    int iters;
    double left;
    double right;
    double lower;
    double upper;
    int width;
    int height;
    double x_part;
    double y_part;
    int *image;
} ARG;

void* cal(void* arg){
    ARG *data = (ARG*)arg;
    int iters = data->iters;
    double left = data->left;
    double lower = data->lower;
    int width = data->width;
    int height = data->height;
    double x_part = data->x_part;
    double y_part = data->y_part;
    int *image = data->image;
    int j ;
    while(1){
        int vec;
        pthread_mutex_lock(&mutex);
        j = j_cur;
        if(j >= height){
            pthread_mutex_unlock(&mutex);
            break;
        }
        
        if(j_cur + 3 < height){
            j_cur += 3;
            vec = 3;
        }
        else if(j_cur + 2 < height){
            j_cur += 2;
            vec = 2;
        }
        else {
            j_cur++;
            vec = 1;
        }
        pthread_mutex_unlock(&mutex);
        
        double y0 = j * y_part + lower;
        double yy0 = (j+1) * y_part + lower;
        double yyy0 = (j+2) * y_part + lower;
        for (int i = 0; i < width; ++i) {
            double x0 = i * x_part + left;
            int repeats = 0;
            double x = 0;
            double y = 0;
            double length_squared = 0;

            double x0_2 = i * x_part + left;
            int repeats_2 = 0;
            double xx = 0;
            double yy = 0;
            double length_squared_2 = 0;
            
            
            double x0_3 = i * x_part + left;
            int repeats_3 = 0;
            double xxx = 0;
            double yyy = 0;
            double length_squared_3 = 0;
            

            while ((repeats < iters && length_squared < 4) || (vec > 1 && repeats_2 < iters && length_squared_2 < 4)|| (vec > 2 && repeats_3 < iters && length_squared_3 < 4)) {
                if(repeats < iters && length_squared < 4){
                    double temp = x * x - y * y + x0;
                    y = 2 * x * y + y0;
                    x = temp;
                    length_squared = x * x + y * y;
                    ++repeats;
                }
                if(vec > 1 && repeats_2 < iters && length_squared_2 < 4){
                    double temp_2 = xx * xx - yy * yy + x0_2;
                    yy = 2 * xx * yy + yy0;
                    xx = temp_2;
                    length_squared_2 = xx * xx + yy * yy;
                    ++repeats_2;
                }
                
                if(vec > 2 && repeats_3 < iters && length_squared_3 < 4){
                    double temp_3 = xxx * xxx - yyy * yyy + x0_3;
                    yyy = 2 * xxx * yyy + yyy0;
                    xxx = temp_3;
                    length_squared_3 = xxx * xxx + yyy * yyy;
                    ++repeats_3;
                }
                
            }
            image[j * width + i] = repeats;
            if(vec > 1) image[(j+1) * width + i] = repeats_2;
            if(vec > 2) image[(j+2) * width + i] = repeats_3;
        }
    }
    pthread_exit(NULL);
}


int main(int argc, char** argv) {
    /* detect how many CPUs are available */
    cpu_set_t cpu_set;
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
    int ncpus = CPU_COUNT(&cpu_set);
    pthread_t thread[ncpus];
    pthread_mutex_init (&mutex, NULL);
    
    /* argument parsing */
    assert(argc == 9);

    ARG arg;
    arg.filename = argv[1];
    arg.iters = strtol(argv[2], 0, 10);
    arg.left = strtod(argv[3], 0);
    arg.right = strtod(argv[4], 0);
    arg.lower = strtod(argv[5], 0);
    arg.upper = strtod(argv[6], 0);
    arg.width = strtol(argv[7], 0, 10);
    arg.height = strtol(argv[8], 0, 10);
    arg.x_part = ((arg.right - arg.left) / arg.width);
    arg.y_part = ((arg.upper - arg.lower) / arg.height);
    /* allocate memory for image */
    int* image = (int*)malloc(arg.width * arg.height * sizeof(int));
    assert(image);
    arg.image = image;

    for(int i = 0; i < ncpus ; i++) pthread_create(&thread[i],NULL,cal,&arg);
    
    void *s;
    for(int i = 0;i<ncpus;i++)
      pthread_join(thread[i],&s);
    
    /* draw and cleanup */
    write_png(arg.filename, arg.iters, arg.width, arg.height, arg.image);
    free(arg.image);
    pthread_mutex_destroy(&mutex);
   	pthread_exit(NULL);
}

