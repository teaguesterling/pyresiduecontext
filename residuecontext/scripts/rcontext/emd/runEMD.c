#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "emd.h"

float norm2 (feature_t *u, feature_t *v) {
  float dx = u->x - v->x,
        dy = u->y - v->y,
        dz = u->z - v->z;
  return sqrt(dx*dx + dy*dy + dz*dz);
}

void usage(const char *arg) {
  printf("Usage %s: sig1 sig2 size\n", arg);
}

void debug(signature_t sig) {
  fprintf(stderr, "SIZE: %d\n", sig.n);
  for(int i = 0; i < sig.n; i++) {
    fprintf(stderr, "%f	%f	%f\n", sig.Features[i].x, sig.Features[i].y, sig.Features[i].z);
  }
  for(int i = 0; i < sig.n; i++) {
    fprintf(stderr, "%f\n", sig.Weights[i]);
  }
}

signature_t load_signature(const char *path, int size) {
  float *weights = (float*) calloc(sizeof(float), size);
  feature_t *features = (feature_t*) calloc(sizeof(feature_t), size);
  FILE *source = fopen(path, "r");
  for(int i = 0; i < size; i++) {
    if(fscanf(source, "%f,%f,%f,%f,", 
                &features[i].x, 
                &features[i].y, 
                &features[i].z, 
                &weights[i]) != 4) {
      fprintf(stderr, "Incomplete signature at item %d\n", i);
    }
  }
  fclose(source);
  signature_t signature = { size, features, weights }; 
  return signature;
}

int main (int argc, const char *argv[]) {
  if(argc != 4) {
    usage(argv[0]);
    return -1;
  }

  int size = atoi(argv[3]);

  signature_t signature1 = load_signature(argv[1], size),
              signature2 = load_signature(argv[2], size);

  #if DEBUG > 1
  debug(signature1);
  debug(signature2);
  #endif

  float dist = emd(&signature1, &signature2, &norm2, 0, 0);

  printf("%f\n", dist);

  return 0;
}
