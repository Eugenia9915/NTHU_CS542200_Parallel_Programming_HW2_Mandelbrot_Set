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

#include <mpi.h>
#include <omp.h>
#include <time.h>

omp_lock_t l;

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

int main(int argc, char** argv) {
    /* argument parsing */
    assert(argc == 9);
    const char* filename = argv[1];
    int iters = strtol(argv[2], 0, 10);
    double left = strtod(argv[3], 0);
    double right = strtod(argv[4], 0);
    double lower = strtod(argv[5], 0);
    double upper = strtod(argv[6], 0);
    int width = strtol(argv[7], 0, 10);
    int height = strtol(argv[8], 0, 10);
    double x_part = ((right - left) / width);
    double y_part = ((upper - lower) / height);
    
    /* detect how many CPUs are available */
    cpu_set_t cpu_set;
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
    int ncpus = CPU_COUNT(&cpu_set);
    //printf("%d cpus available\n", ncpus);
    
    int rank, nproc;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);

    int remain = height % nproc, j_size_per_p;
    int i = 0 , j = 0, j_start , j_end;
    if(rank < remain){
        j_size_per_p = height / nproc + 1;
        j = rank * j_size_per_p;
        j_start = j;
        j_end = j + j_size_per_p;
    }
    else{
        j_size_per_p = height / nproc;
        j = remain + rank * j_size_per_p;
        j_start = j;
        j_end = j + j_size_per_p;
    }

    omp_init_lock(&l);

    /* allocate memory for image */
    int* image = (int*)malloc(width * j_size_per_p * sizeof(int));
    assert(image);
    int num = 0;
    #pragma omp parallel num_threads(ncpus)
    {
        int ii, jj, n;
        int vec;
        while(1){
            omp_set_lock(&l);
            jj = j; 
            n = num;
            if(jj == j_end){
                omp_unset_lock(&l);
                break;
            }
            if(j + 3 < j_end){
                j += 3;
                num += 3;
                vec = 3;
            }
            else if(j + 2 < j_end){
                j += 2;
                num += 2;
                vec = 2;
            }
            else {
                j ++;
                num ++;
                vec = 1;
            }
            omp_unset_lock(&l);

            double y0 = jj * y_part + lower;
            double yy0 = (jj+1) * y_part + lower;
            double yyy0 = (jj+2) * y_part + lower;
            for(int ii = 0;ii < width;ii++){
                double x0 = ii * x_part + left;
                int repeats = 0;
                double x = 0;
                double y = 0;
                double length_squared = 0;

                double x0_2 = ii * x_part + left;
                int repeats_2 = 0;
                double xx = 0;
                double yy = 0;
                double length_squared_2 = 0;

                double x0_3 = ii * x_part + left;
                int repeats_3 = 0;
                double xxx = 0;
                double yyy = 0;
                double length_squared_3 = 0;
                while ((repeats < iters && length_squared < 4)||(vec > 1 && repeats_2 < iters && length_squared_2 < 4)|| (vec > 2 && repeats_3 < iters && length_squared_3 < 4)) {
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
                image[(jj-j_start) * width + ii] = repeats;
                if(vec > 1) image[(jj+1-j_start) * width + ii] = repeats_2;
                if(vec > 2) image[(jj+2-j_start) * width + ii] = repeats_3;
            }
        }
    }

    if(rank == 0){
        //處理rank = 0 的 image
        int* image_w = (int*)malloc(width * height * sizeof(int));
        assert(image_w);

        int index = 0;
        for(int i=0;i<width * j_size_per_p;i++,index++){
            image_w[index] = image[i];
        }
        free(image);
        for(int i=1;i<nproc;i++){
            if(i<remain){
                int* image_tmp = (int*)malloc((height / nproc + 1) * width * sizeof(int));
                assert(image_tmp);
                MPI_Recv(image_tmp,(height / nproc + 1)*width,MPI_INT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                for(int j=0;j<(height / nproc + 1)*width;j++,index++)
                image_w[index] = image_tmp[j];
                free(image_tmp);
            }
            else{
                int* image_tmp = (int*)malloc((height / nproc) * width * sizeof(int));
                assert(image_tmp);
                MPI_Recv(image_tmp,(height / nproc)*width,MPI_INT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                for(int j=0;j<(height / nproc)*width;j++,index++)
                image_w[index] = image_tmp[j];
                free(image_tmp);
            }
        }
        /* draw and cleanup */
        write_png(filename, iters, width, height, image_w);
        free(image_w);
    }
    else{
        MPI_Send(image,width * j_size_per_p,MPI_INT,0,1,MPI_COMM_WORLD);
        free(image);
    }
    MPI_Finalize();
}

