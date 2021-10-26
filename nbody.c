#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define N 4 // num children nodes
#define DOMAIN_SIZE 1. // total size of domain
#define THETA 1
#define SOFTENING 1e-9f

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
//void prune(QTNode *node);
void postorder(QTNode *node);
void ran_init(float *data, int n);
void qTree_destroy(QTNode *node);
float *qTree_force(Particle *k, QTNode *n);

int main(const int argc, const char** argv) {
    
    // two tests for QTree design
    
    // Set DOMAIN_SIZE to 16.
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

    // Set DOMAIN_SIZE to 1.
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
    
    FILE *datafile    = NULL;      // output file for particle positions 
    int   nParticles  = 16;      // number of particles 

    if (argc > 1)
        nParticles      = atoi(argv[1]);

    const float dt    = 0.01f; // time step  
    // const int nIters  = 200;   // number of steps in simulation 

    float *buf        =  malloc(nParticles*sizeof(Particle));
    Particle  *p          = (Particle *) buf;
    ran_init(buf, 4*nParticles); // Init pos and vel data 

    QTNode *root = qTree_build(p,nParticles);
    //printf("%f, %f, %f", root->total_mass, root->center_of_mass[0],root->center_of_mass[1]);
    /*for (int i = 0; i < N; i++){
        if (root->child[i] != NULL){
            printf("%d, %f, %f\n", root->child[i]->which_child, root->child[i]->particle->x,root->child[i]->particle->y);
        }
    } */

    datafile          = fopen("particles.dat","w");
    for (int i = 0;i < nParticles; ++i)
        fprintf(datafile, "%f %f \n", p[i].x, p[i].y);
    fclose(datafile);

    for (int i = 0; i < nParticles; i++){
        float *forces_particle = qTree_force(&p[i],root);
        p[i].vx += dt*forces_particle[0];
        p[i].vy += dt*forces_particle[1];
    } 
    
    free(buf);
    qTree_destroy(root); 


    // tests for algorithms accuracy
    
    // Set DOMAIN_SIZE to 16.

    
    return 0;
}

float *qTree_force(Particle *k, QTNode *n){
    static float forces[2] = {0.0f, 0.0f};
    //if (n == NULL) // || n->total_mass == 0)
        //return NULL;
        //printf("N is NULL");
    //float force = 0.;
    //float Fx = 0.0f; float Fy = 0.0f;
    if (n->total_mass == 1){
        float dx = n->particle->x - k->x;
        float dy = n->particle->y - k->y;
        float distSqr = dx*dx + dy*dy + SOFTENING;
        float invDist = 1.0f / sqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;
        float Fx = dx * invDist3; 
        float Fy = dy * invDist3;
        //k->vx += dt*Fx; k->vy += dt*Fy;
        forces[0] += Fx;
        forces[1] += Fy;
        return forces;
    }
    else{
        float r = 0;
        r += sqrtf((k->x - n->center_of_mass[0])*(k->x - n->center_of_mass[0]) + (k->y - n->center_of_mass[1])*(k->y - n->center_of_mass[1])); // r is distance from particle k to CM of particles in n
        if (r != 0 && (n->size / r < THETA)){
            float dx = n->center_of_mass[0] - k->x;
            float dy = n->center_of_mass[1] - k->y;
            float distSqr = dx*dx + dy*dy + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;
            float Fx = n->total_mass * dx * invDist3; 
            float Fy = n->total_mass * dy * invDist3;
            forces[0] += Fx;
            forces[1] += Fy;
            //k->vx += dt*Fx; k->vy += dt*Fy;
            return forces;
        }
        else{ 
            for (int i = 0; i < N; i++){
                if (n->child[i] != NULL){
                    float *force_child = qTree_force(k,n->child[i]);
                    forces[0] += force_child[0];
                    forces[1] += force_child[1];
                }
            }
            return forces;
        }
    }
}

QTNode *qTree_build(Particle *arr, int num_particles){
    QTNode *root = create_node(NULL,0);
    for (int j = 0; j < num_particles; j++){ 
        qTree_insert(&arr[j],root); // arr[j] = pointer to jth particle
    }
    compute_COM(root);
    // prune
    return root;
}

void qTree_insert(Particle *p, QTNode *n){
    QTNode *c;
    /*if (n == NULL)
        printf("n is null"); */
    if (!is_leaf(n)){ // node has 4 children
        // printf("If\n");
        c = which_child_contains(n,p);
        qTree_insert(p,c);
    }
    else if (n->particle != NULL){// node is a leaf that contains a particle
        // printf("Else if\n");
        for (int i=0; i < N; i++){
            // printf("Children created\n");
            n->child[i] = create_node(n,i);
        }
        // printf("%f,%f",n->particle->x,n->particle->y);
        c = which_child_contains(n,n->particle); // find the child which contains particle p in the parent node
        /* if (c == NULL)
            printf("c is NULL");
        else
            printf("%d", c->which_child); */
        c->particle = n->particle;
        c->total_mass++;
        c->center_of_mass[0] = c->particle->x;
        c->center_of_mass[1] = c->particle->y;
        //printf("%f",c->particle->x);
        n->particle = NULL;
        n->total_mass--;
        n->center_of_mass[0] = 0; // initialized to address warning from Valgrind in qtree_force
        n->center_of_mass[1] = 0; // initialized to address warning from Valgrind in qtree_force
        c = which_child_contains(n,p);
        /*if (c == NULL)
            printf("c is NULL");
        else
            printf("%d", c->which_child); */
        qTree_insert(p,c);
    }
    else{ // node is a leaf that contains no particles
        // printf("Else, %f %f\n",p->x, p->y);
        n->particle = p;
        n->total_mass++;

        n->center_of_mass[0] = n->particle->x;
        n->center_of_mass[1] = n->particle->y;
        /*QTNode *node = n->parent;
        while(node != NULL) {
            //n->parent->center_of_mass[0] = ((n->parent->center_of_mass[0]) * (n->parent->total_mass) + (n->particle->x)) / (n->parent->total_mass + 1);
            //n->parent->center_of_mass[1] = ((n->parent->center_of_mass[1]) * (n->parent->total_mass) + (n->particle->y)) / (n->parent->total_mass + 1);
            node->total_mass++;
            node = node->parent;
        } */
    }
}

void compute_COM(QTNode *root){
    if (root == NULL || is_leaf(root)) 
        return;
    else {
        for (int i = 0; i < N; i++)
            compute_COM(root->child[i]);
        float x_coord_total = 0;
        float y_coord_total = 0;
        for (int i = 0; i < N; i++){
            root->total_mass += root->child[i]->total_mass;
            x_coord_total += root->child[i]->center_of_mass[0] * root->child[i]->total_mass;
            y_coord_total += root->child[i]->center_of_mass[1] * root->child[i]->total_mass;
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
    node->center_of_mass[0] = 0; // legal?
    node->center_of_mass[1] = 0;
    if (parent == NULL){
        node->size = DOMAIN_SIZE;
        node->lb = 0.; node->rb = DOMAIN_SIZE;
        node->ub = DOMAIN_SIZE; node->db = 0.; 
        /* node->size = 16.;
        node->lb = 0.; node->rb = 16.;
        node->ub = 16.; node->db = 0.; */
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
            node->rb = node->parent->lb + ((node->parent->rb - node->parent->lb)/float_two);//node->parent->rb/float_two;
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
            node->rb = node->parent->lb + ((node->parent->rb - node->parent->lb)/float_two);//node->parent->rb/float_two;
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
        //if (p->x < n->child[i]->rb && p->x > n->child[i]->lb && p->y < n->child[i]->ub && p->y > n->child[i]->db)
        //    return n->child[i];
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

void postorder(QTNode *node){
    if (node == NULL)
        return;
    for (int i = 0; i <= 3; i++)
        postorder(node->child[i]);
    // if (node->particle != NULL)
    //    printf("[%f,%f]\n",node->p..) */
}

/* randomly initialize particle positions and momenta */
void ran_init(float *data, int n) {
  for (int i = 0; i < n; i++) {
    data[i] = (float)rand()/(float)(RAND_MAX);
    // data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
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

/*void qTree_print(){
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            printf("%.2f,%.2f ",root->child[i]->child[j]->particle->x,root->child[i]->child[j]->particle->y);
        }
        printf("\n");
    } 
} */

/* void calc_force(Particle *p, float dt, int n) {
  for (int i = 0; i < n; i++) { 
    float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

    for (int j = 0; j < n; j++) {
      float dx = p[j].x - p[i].x;
      float dy = p[j].y - p[i].y;
      float dz = p[j].z - p[i].z;
      float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
    }
    p[i].vx += dt*Fx; p[i].vy += dt*Fy; p[i].vz += dt*Fz;
  }
} */
