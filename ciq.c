// K-means++ clustering algorithm for 8-bit color images
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Uncomment the following line to enable debug mode
// #define __DEBUG__
#define MAX_ITERS 100   // maximum number of iterations
#define EPSILON 8       // threshold for centroid update

// KCIQ: Define boolean type
#ifndef bool
    #define bool int
    #define true 1
    #define false 0
#endif

// KCIQ: Define structure for data points and centroids
typedef struct {
    int r, g, b;
    int cluster;
} Point;

typedef Point Centroid;

typedef struct context {
    int width, height;
    long size;
    int K;
    Point * points;
    Centroid * centroids;
} Context;

// KCIQ: calculate Euclidean distance
long ciq_distance(Point p1, Centroid p2) {
    long dr = p1.r - p2.r;
    long dg = p1.g - p2.g;
    long db = p1.b - p2.b;
    return (dr*dr + dg*dg + db*db);
}

// KCIQ: initialize the context
Context * ciq_init(const char * filename, int K) {
    Context * ctx = (Context *) malloc(sizeof(Context));
    if (!ctx) {
#ifdef  __DEBUG__
        fprintf(stderr, "Memory allocation failed\n");
#endif
        return false;
    }

    // Open the file
    FILE * file = fopen(filename, "rb");
    if (!file) {
#ifdef __DEBUG__
        fprintf(stderr, "Unable to open file %s\n", filename);
#endif
        free(ctx);
        return NULL;
    }

    // Read PPM header, please note that we does not support comments here
    char format[3];
    int width, height, maxval;
    fscanf(file, "%s\n%d %d\n%d\n", format, &width, &height, &maxval);

    if (strcmp(format, "P6") != 0 || maxval != 255) {
#ifdef __DEBUG__        
        fprintf(stderr, "Unsupported PPM format\n");
#endif
        free(ctx);
        fclose(file);
        return false;
    }

    // update the context
    ctx->width = width;
    ctx->height = height;
    ctx->size = width * height;
    ctx->K = K;

#ifdef __DEBUG__
    printf("- Image size: %dx%d\n", width, height);
    printf("- Number of data points: %ld\n", ctx->size);
    printf("- Number of clusters: %d\n", K);
#endif

    // allocate memory for data points and centroids
    ctx->points = (Point *) malloc(ctx->size * sizeof(Point));
    ctx->centroids = (Centroid *) malloc(K * sizeof(Centroid));
    if (!ctx->points || !ctx->centroids) {
#ifdef __DEBUG__
        fprintf(stderr, "Memory allocation failed\n");
#endif
        if (ctx->points) free(ctx->points);
        if (ctx->centroids) free(ctx->centroids);
        free(ctx);
        fclose(file);
        return NULL;
    }
    else {
        memset(ctx->points, 0, ctx->size * sizeof(Point));
        memset(ctx->centroids, 0, K * sizeof(Centroid));
#ifdef  __DEBUG__
        printf("- Allocated %lu bytes for the data points\n", ctx->size * sizeof(Point));
        printf("- Allocated %lu bytes for the centroids\n", K * sizeof(Centroid));
#endif        
    }

    // read the image data
    for (int i = 0; i < ctx->size; i++) {
        unsigned char r, g, b;
        fread(&r, 1, 1, file);
        fread(&g, 1, 1, file);
        fread(&b, 1, 1, file);
        ctx->points[i] = (Point) {r, g, b, -1};
    }

    fclose(file);   // close the file
    return ctx;     // return the context
}

// KCIQ: K-means++ initialization
bool ciq_init_centroids(Context * ctx) {

    if (!ctx) return false;

    int i, j;
    int chosen_index;
    long total_distance, random_choice, cumulative_probability;
    long *distances = (long *) malloc(ctx->size * sizeof(long));

    if (!distances)
        return false;

    // Choose the first centroid randomly
    chosen_index = rand() % ctx->size;
    ctx->centroids[0] = (Centroid){ ctx->points[chosen_index].r, 
                                    ctx->points[chosen_index].g, 
                                    ctx->points[chosen_index].b};
#ifdef __DEBUG__
    printf("- Initial centroid: (%d, %d, %d)\n", 
            ctx->centroids[0].r, ctx->centroids[0].g, ctx->centroids[0].b);
#endif

    // Choose the remaining centroids
    for (i = 1; i < ctx->K; i++) {
        total_distance = 0.0;
        for (j = 0; j < ctx->size; j++) {
            distances[j] = ciq_distance(ctx->points[j], ctx->centroids[i - 1]);
            total_distance += distances[j];
        }

        random_choice = ((double) rand() / RAND_MAX) * total_distance;
        cumulative_probability = 0.0;
        for (j = 0; j < ctx->size; j++) {
            cumulative_probability += distances[j];
            if (cumulative_probability >= random_choice) {
                ctx->centroids[i] = (Centroid){ ctx->points[j].r, 
                                                ctx->points[j].g, 
                                                ctx->points[j].b};
#ifdef __DEBUG__
                printf("- Centroid %3d: (%d, %d, %d)\n", i,
                        ctx->centroids[i].r, ctx->centroids[i].g, ctx->centroids[i].b);
#endif                                                
                break;
            }
        }
    }
    free(distances);
    return true;
}

// KCIQ: assign points to the nearest centroid
void ciq_clustering(Context * ctx) {
    int i, j;
    long mindist, curdist;

    if (!ctx) return;    
    for (i = 0; i < ctx->size; i++) {
        mindist = ciq_distance(ctx->points[i], ctx->centroids[0]);
        ctx->points[i].cluster = 0;
        for (j = 1; j < ctx->K; j++) {
            curdist = ciq_distance(ctx->points[i], ctx->centroids[j]);
            if (curdist < mindist) {
                mindist = curdist;
                ctx->points[i].cluster = j;
            }
        }
    }
}

// KCIQ: update centroids based on assigned points
bool ciq_update_centroids(Context * ctx) {
    
    if (!ctx) return false;
    
    Centroid new;
    int i, j, cluster_size[ctx->K];
    long sum_r[ctx->K], sum_g[ctx->K], sum_b[ctx->K];
    double w[ctx->K];
    bool changed = false;

    // reset the sums and cluster sizes
    memset(sum_r, 0, sizeof(sum_r));
    memset(sum_g, 0, sizeof(sum_g));
    memset(sum_b, 0, sizeof(sum_b));
    memset(cluster_size, 0, sizeof(cluster_size));

    // calculate the sums and cluster sizes
    for (i = 0; i < ctx->size; i++) {
        j = ctx->points[i].cluster;
        sum_r[j] += ctx->points[i].r;
        sum_g[j] += ctx->points[i].g;
        sum_b[j] += ctx->points[i].b;
        cluster_size[j]++;
    }

    // calculate the factor per cluster
    for (i = 0; i < ctx->K; i++)
        w[i] = 1.0 / cluster_size[i];

    // update the centroids
    for (i = 0; i < ctx->K; i++) {
        if (cluster_size[i] > 0) {
            new.r = w[i] * sum_r[i];
            new.g = w[i] * sum_g[i];
            new.b = w[i] * sum_b[i];
        }
        // check if the centroid has changed
        if (ciq_distance(ctx->centroids[i], new) > EPSILON) {
            changed = true;
        }
        // update the current centroid
        ctx->centroids[i] = new;
    }
    return changed;
}

// KCIQ: free memory
void ciq_shutdown(Context * ctx) {
    if (!ctx) return;
    if (ctx->points) 
        free(ctx->points);
    if (ctx->centroids)
        free(ctx->centroids);
    free(ctx);
}

// KCIQ: perform k-means clustering for image quantization
bool ciq_quantize(Context * ctx) {
    int i;

    if (!ctx) return false;    
    if (!ciq_init_centroids(ctx))
        return false;

    for (i = 0; i < MAX_ITERS; i++) {  
        printf("Iteration: %d\r", i+1);
        ciq_clustering(ctx);
        bool changed = ciq_update_centroids(ctx);
        if (!changed) {
#ifdef __DEBUG__
        if (!changed)
            printf("\n- Clusters stable.\n");
#endif            
            break;
        }
        fflush(stdout);
    }
    printf("\n");
    return true;
}

// KCIQ: remap the image using the quantized palette
bool ciq_remap(Context * ctx, const char * filename) {
    if (!ctx) return false;

    FILE *file = fopen(filename, "wb");
    if(!file) {
#ifdef __DEBUG__
        fprintf(stderr, "Unable to create file %s\n", filename);
#endif
        return false;
    }

    fprintf(file, "P6\n%d %d\n255\n", ctx->width, ctx->height);

    for (int i = 0; i < ctx->size; i++) {
        Centroid c = ctx->centroids[ctx->points[i].cluster];
        unsigned char r = c.r;
        unsigned char g = c.g;
        unsigned char b = c.b;
        fwrite(&r, 1, 1, file);
        fwrite(&g, 1, 1, file);
        fwrite(&b, 1, 1, file);
    }

    fclose(file);

    // write the palette file
    file = fopen("palette.pal", "wb");
    if (!file) {
#ifdef __DEBUG__
        fprintf(stderr, "Unable to create file palette.pal\n");
#endif
        return false;
    }

    for (int i = 0; i < ctx->K; i++) {
        Centroid c = ctx->centroids[i];
        unsigned char r = c.r;
        unsigned char g = c.g;
        unsigned char b = c.b;
        fwrite(&r, 1, 1, file);
        fwrite(&g, 1, 1, file);
        fwrite(&b, 1, 1, file);
    }
    fclose(file);
    return true;
}

// KCIQ: main function for image quantization
bool ciq_quanization(const char * input, const char * output, int K) {
    Context * ctx = ciq_init(input, K);
    if (!ctx) {
#ifdef __DEBUG__
        fprintf(stderr, "Failed to initialize context\n");
#endif        
        return false;
    }
    ciq_quantize(ctx);
    if (ciq_remap(ctx, output)) {
#ifdef __DEBUG__
        printf("Image quantized successfully and saved into %s\n", output);
#endif        
    } else {
#ifdef  __DEBUG__        
        fprintf(stderr, "Failed to quantize image\n");
#endif        
    }
    ciq_shutdown(ctx);
    return true;
}

int main(int argc, char *argv[]) {
    printf("Color Image Quantization using K-Means++ - v0.1\n");
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input.ppm> <output.ppm> [K]\n", argv[0]);
        return 1;
    }

    char input[256];
    char output[256];
    int K = argv[3] ? atoi(argv[3]) : 256;
    
    strcpy(input, argv[1]);
    strcpy(output, argv[2]);

    printf("Quantizing image %s with K=%d\n", input, K);

    if (!ciq_quanization(input, output, K)) {
        fprintf(stderr, "Failed to quantize image\n");
        return 1;
    }

    return 0;
}