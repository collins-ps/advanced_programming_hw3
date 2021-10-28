#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "timer.h"

#define N 4 // num children nodes
#define DOMAIN_SIZE 1. // total size of domain
#define THETA 0
#define SOFTENING 1e-9f
#define MARGIN 0.01

typedef struct particle_{
  float x, y;        /* particle positions */
  float vx, vy;     /* particle momenta */
} Particle;

typedef struct qtnode_{
    Particle *particle;
    int which_child;
    float size; // size of domain represented by the node
    float center_of_mass[2]; // spatial average position (assuming unit mass)
    float total_mass; // number of bodies, assuming unit mass
    float lb, rb, db, ub; // physical space domain boundary
    struct qtnode_* parent;
    struct qtnode_* child[N]; // pointer to N children nodes
} QTNode;

QTNode *qTree_build(Particle *arr, int num_particles);
void qTree_insert(Particle *p, QTNode *n);
void compute_COM(QTNode *root);
QTNode *create_node(QTNode *parent, int child_index);
QTNode *which_child_contains(QTNode *n, Particle *p);
static int in_nox(Particle *p, float lb, float rb, float db, float ub);
int is_leaf(QTNode *n);
void ran_init(float *data, int n);
void qTree_destroy(QTNode *node);
float *qTree_force(Particle *k, QTNode *n);
void calc_force(Particle *p, float dt, int n);

int main(const int argc, const char** argv) {
    
    FILE *datafile    = NULL;      // output file for particle positions 
    int   nParticles  = 3000;      // number of particles 

    if (argc > 1)
        nParticles      = atoi(argv[1]);

    const float dt    = 0.01f; // time step  
    const int nIters  = 1;   // number of steps in simulation 

    float *buf        =  malloc(nParticles*sizeof(Particle));
    Particle  *p          = (Particle *) buf;
    ran_init(buf, 4*nParticles); // Init pos and vel data 
    
    // create copy of buf for testing purposes //
    float *buf_test        =  malloc(nParticles*sizeof(Particle));
    Particle  *p_test          = (Particle *) buf_test;
    for (int k = 0; k < 4*nParticles; k++){
        buf_test[k] = buf[k];
    }

    double totalTime  = 0.0;

    datafile          = fopen("particles.dat","w");
    fprintf(datafile,"%d %d %d\n", nParticles, nIters, 0);

    /* ------------------------------*/
    /*     MAIN LOOP                 */
    /* ------------------------------*/
    QTNode *root = qTree_build(p,nParticles);

    for (int iter = 1; iter <= nIters; iter++) {
        printf("iteration:%d\n", iter);
    
        for (int i = 0;i < nParticles; ++i)
            fprintf(datafile, "%f %f \n", p[i].x, p[i].y);
        
        StartTimer();

        for (int i = 0; i < nParticles; i++){
            float *forces_particle = qTree_force(&p[i],root);
            p[i].vx += dt*forces_particle[0];
            p[i].vy += dt*forces_particle[1];
            free(forces_particle);
        } 
        for (int i = 0 ; i < nParticles; i++) {  /* compute new position */
            p[i].x += p[i].vx*dt;
            p[i].y += p[i].vy*dt;
        }
        
        const double tElapsed = GetTimer() / 1000.0;
        if (iter > 1) {                          /* First iter is warm up */
            totalTime += tElapsed; 
        }

        // tests for algorithm's accuracy
        calc_force(p_test, dt, nParticles); 
        for (int i = 0 ; i < nParticles; i++) {  /* compute new position */
            p_test[i].x += p_test[i].vx*dt;
            p_test[i].y += p_test[i].vy*dt;
        }
        /*for (int i = 0 ; i < nParticles; i++) { 
            if((fabs(p[i].x - p_test[i].x) > MARGIN) || (fabs(p[i].y - p_test[i].y) > MARGIN))
                printf("Error: x and y differ by %f and %f for particle %d.\n", fabs(p[i].x - p_test[i].x),fabs(p[i].y - p_test[i].y),i); 
        } */
        for (int i = 0 ; i < nParticles; i++) { 
            assert(fabs(p[i].x - p_test[i].x) < MARGIN); 
            assert(fabs(p[i].y - p_test[i].y) < MARGIN); 
        }
        printf("Passed accuracy tests.\n"); 
        
    }

    fclose(datafile);
    double avgTime = totalTime / (double)(nIters-1); 
    printf("avgTime: %f   totTime: %f \n", avgTime, totalTime);

    qTree_destroy(root); 
    free(buf);
    free(buf_test);

    return 0;

    // two tests for QTree design
    
    // #define DOMAIN_SIZE 16.
    /*
    int nParticles  = 16;      
    float buf[64] = {2,14,0,0,6,14,0,0,9,14,0,0,14,14,0,0,2,9,0,0,6,9,0,0,9,9,0,0,14,9,0,0,2,6,0,0,6,6,0,0,9,6,0,0,14,6,0,0,2,2,0,0,6,2,0,0,9,2,0,0,14,2,0,0};
    Particle *p = (Particle *) buf;
    QTNode *root = qTree_build(p,nParticles);
    for (int i = 0; i < 4; i++){
        assert(root->child[i]->total_mass == 4);
        for (int j = 0; j < 4; j++){
            assert(root->child[i]->child[j]->particle != NULL);
            assert(root->child[i]->child[j]->which_child == j);
            assert(root->child[i]->child[j]->size == DOMAIN_SIZE/16);
            assert(root->child[i]->child[j]->total_mass == 1);
            assert(root->child[i]->child[j]->parent == root->child[i]);
            assert(root->child[i]->child[j]->center_of_mass[0] == root->child[i]->child[j]->particle->x);
            assert(root->child[i]->child[j]->center_of_mass[1] == root->child[i]->child[j]->particle->y);
        }
    }
    assert(root->size == 16.);
    assert(root->total_mass == 16);
    assert(root->center_of_mass[0] == 7.75);
    assert(root->center_of_mass[1] == 7.75);
    assert(root->child[0]->center_of_mass[0] == 4);
    assert(root->child[0]->center_of_mass[1] == 11.5);
    assert(root->child[1]->center_of_mass[0] == 11.5);
    assert(root->child[1]->center_of_mass[1] == 11.5);
    assert(root->child[2]->center_of_mass[0] == 4);
    assert(root->child[2]->center_of_mass[1] == 4);
    assert(root->child[3]->center_of_mass[0] == 11.5);
    assert(root->child[3]->center_of_mass[1] == 4);

    printf("Passed all tests.\n");
    qTree_destroy(root); 
    */

    // #define DOMAIN_SIZE 1.
    /*
    int nParticles  = 16;      
    float buf[64] = {0.2,0.8,0,0,0.3,0.8,0,0,0.6,0.8,0,0,0.8,0.8,0,0,0.2,0.6,0,0,0.3,0.6,0,0,0.6,0.6,0,0,0.8,0.6,0,0,0.2,0.3,0,0,0.3,0.3,0,0,0.6,0.3,0,0,0.8,0.3,0,0,0.2,0.2,0,0,0.3,0.2,0,0,0.6,0.2,0,0,0.8,0.2,0,0};
    Particle *p = (Particle *) buf;
    QTNode *root = qTree_build(p,nParticles);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            assert(root->child[i]->child[j]->particle != NULL);
            assert(root->child[i]->child[j]->which_child == j);
            assert(root->child[i]->child[j]->size == DOMAIN_SIZE/16);
            assert(root->child[i]->child[j]->total_mass == 1);
            assert(root->child[i]->child[j]->parent == root->child[i]);
        }
    }
    printf("Passed all tests.\n");
    qTree_destroy(root); 
    */
}

float *qTree_force(Particle *k, QTNode *n){
    float *forces = malloc(2*sizeof(float));
    forces[0] = 0; forces[1] = 0;
    if (n->total_mass == 1){
        float dx = n->particle->x - k->x;
        float dy = n->particle->y - k->y;
        float distSqr = dx*dx + dy*dy + SOFTENING;
        float invDist = 1.0f / sqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;
        forces[0] = dx * invDist3; 
        forces[1] = dy * invDist3;
        return forces;
    }
    else{
        float r = sqrtf((k->x - n->center_of_mass[0])*(k->x - n->center_of_mass[0]) + (k->y - n->center_of_mass[1])*(k->y - n->center_of_mass[1])); // r is distance from particle k to CM of particles in n
        if (r != 0 && (n->size / r < THETA)){
            float dx = n->center_of_mass[0] - k->x;
            float dy = n->center_of_mass[1] - k->y;
            float distSqr = dx*dx + dy*dy + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;
            forces[0] = n->total_mass * dx * invDist3; 
            forces[1] = n->total_mass * dy * invDist3;
            return forces;
        }
        else{ 
            for (int i = 0; i < N; i++){
                if (n->child[i] != NULL){
                    float *force_child = qTree_force(k,n->child[i]);
                    forces[0] += force_child[0];
                    forces[1] += force_child[1];
                    free(force_child);
                }
            }
            return forces;
        }
    }
}

QTNode *qTree_build(Particle *arr, int num_particles){
    QTNode *root = create_node(NULL,0);
    for (int j = 0; j < num_particles; j++){ 
        qTree_insert(&arr[j],root); 
    }
    compute_COM(root);
    return root;
}

void qTree_insert(Particle *p, QTNode *n){
    QTNode *c;
    if (!is_leaf(n)){ // node has 4 children
        c = which_child_contains(n,p);
        qTree_insert(p,c);
    }
    else if (n->particle != NULL){// node is a leaf that contains a particle
        for (int i=0; i < N; i++){
            n->child[i] = create_node(n,i);
        }
        c = which_child_contains(n,n->particle); // find the child which contains particle p in the parent node
        c->particle = n->particle;
        c->total_mass++;
        n->particle = NULL;
        n->total_mass--;
        c = which_child_contains(n,p);
        qTree_insert(p,c);
    }
    else{ // node is a leaf that contains no particles
        n->particle = p;
        n->total_mass++;
    }
}

void compute_COM(QTNode *root){
    if (root == NULL) 
        return;
    else if (root->particle != NULL && is_leaf(root)){
        root->center_of_mass[0] = root->particle->x;
        root->center_of_mass[1] = root->particle->y;
        return;
    }
    else {
        for (int i = 0; i < N; i++)
            compute_COM(root->child[i]);
        float x_coord_total = 0;
        float y_coord_total = 0;
        for (int i = 0; i < N; i++){
            if (root->child[i] != NULL){
                root->total_mass += root->child[i]->total_mass;
                x_coord_total += root->child[i]->center_of_mass[0] * root->child[i]->total_mass;
                y_coord_total += root->child[i]->center_of_mass[1] * root->child[i]->total_mass;
            }
        }
        root->center_of_mass[0] = x_coord_total / root->total_mass;
        root->center_of_mass[1] = y_coord_total / root->total_mass;
        return;
    }
}

QTNode *create_node(QTNode *parent, int child_index){
    QTNode *node = malloc(sizeof(QTNode));
    node->which_child = child_index;
    for (int i = 0; i < N; i++){
        node->child[i] = NULL;
    }
    node->particle = NULL;
    node->parent = parent;
    node->total_mass = 0;
    if (parent == NULL){
        node->size = DOMAIN_SIZE;
        node->lb = 0.; node->rb = DOMAIN_SIZE;
        node->ub = DOMAIN_SIZE; node->db = 0.; 
    }
    else{
        float float_four = 4.0;
        float float_two = 2.0;
        node->size = parent->size/float_four;
        if (node->which_child == 0){
            // y coordinates
            node->ub = node->parent->ub;
            node->db = node->parent->db + ((node->parent->ub - node->parent->db)/float_two);
            // x coordinates
            node->lb = node->parent->lb;
            node->rb = node->parent->lb + ((node->parent->rb - node->parent->lb)/float_two);
        }
        else if (node->which_child == 1){
            // y coordinates
            node->ub = node->parent->ub;
            node->db = node->parent->db + ((node->parent->ub - node->parent->db)/float_two);
            // x coordinates
            node->lb = node->parent->lb + ((node->parent->rb - node->parent->lb)/float_two);
            node->rb = node->parent->rb;
        }
        else if (node->which_child == 2){
            // y coordinates
            node->ub = node->parent->db + ((node->parent->ub - node->parent->db)/float_two);
            node->db = node->parent->db;
            // x coordinates
            node->lb = node->parent->lb;
            node->rb = node->parent->lb + ((node->parent->rb - node->parent->lb)/float_two);
        }
        else if (node->which_child == 3){
            // y coordinates
            node->ub = node->parent->db + ((node->parent->ub - node->parent->db)/float_two);
            node->db = node->parent->db;
            // x coordinates
            node->lb = node->parent->lb + ((node->parent->rb - node->parent->lb)/float_two);
            node->rb = node->parent->rb;
        }
    }
    return node;
}

QTNode* which_child_contains(QTNode *n, Particle *p){
    for (int i = 0; i < N; i++){
        if (in_nox(p,n->child[i]->lb,n->child[i]->rb,n->child[i]->db,n->child[i]->ub))
            return (n->child[i]);
    }
    printf("Error - particle not in node's domain");
    return NULL;
}

static int in_nox(Particle *p, float lb, float rb, float db, float ub){
    if (p->x >= lb && p->x <= rb)
        if (p->y >= db && p->y <= ub)
            return 1;
    return 0;
}

int is_leaf(QTNode *n){
    if (n->child[0] == NULL && n->child[1] == NULL && n->child[2] == NULL && n->child[3] == NULL)
        return(1);
    else
        return(0);
}

/* randomly initialize particle positions and momenta */
void ran_init(float *data, int n) {
  for (int i = 0; i < n; i++) {
    data[i] = fabs(2.0f * (rand() / (float)RAND_MAX) - 1.0f);
    //data[i] = (float)rand()/(float)(RAND_MAX);
  }
}

void qTree_destroy(QTNode *node){
    if (node == NULL) 
        return;
    for (int i = 0; i < 4; ++i) {
        qTree_destroy(node->child[i]);
        node->child[i] = NULL;
    }
    free(node);
}

void calc_force(Particle *p, float dt, int n) {
  for (int i = 0; i < n; i++) { 
    float Fx = 0.0f; float Fy = 0.0f;

    for (int j = 0; j < n; j++) {
      float dx = p[j].x - p[i].x;
      float dy = p[j].y - p[i].y;
      float distSqr = dx*dx + dy*dy + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3; Fy += dy * invDist3; 
    }
    p[i].vx += dt*Fx; p[i].vy += dt*Fy; 
  }
} 
