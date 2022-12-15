#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"


int n_body;
int n_iteration;

int n_omp_threads;


void generate_data(double *m, double *x,double *y,double *vx,double *vy, int n) {
    // TODO: Generate proper initial position and mass for better visualization
    srand((unsigned)time(NULL));
    for (int i = 0; i < n; i++) {
        m[i] = rand() % max_mass + 1.0f;
        x[i] = 2000.0f + rand() % (bound_x / 4);
        y[i] = 2000.0f + rand() % (bound_y / 4);
        vx[i] = 0.0f;
        vy[i] = 0.0f;
    }
}



void update_position(double *x, double *y, double *vx, double *vy, int i) {
    //TODO: update position
    x[i] += vx[i] * dt;
    y[i] += vy[i] * dt;

    if (x[i] > 2000.0f + bound_x / 4)
    x[i] = 2000.0f + bound_x / 4 - 1;
    if (y[i] > 2000.0f + bound_y / 4)
    y[i] = 2000.0f + bound_y / 4 - 1;          
    if (x[i] < 2000.0f - bound_x / 4)
    x[i] = 2000.0f - bound_x / 4 + 1;
    if (y[i] < 2000.0f - bound_y / 4)
    y[i] = 2000.0f - bound_y / 4 + 1;   
}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int pi, int n) {
    //TODO: calculate force and acceleration, update velocity
    double ori_vx;
    double ori_vy;
    ori_vx = vx[pi];
    ori_vy = vy[pi];

    for (int pj = 0; pj < n; pj++){
        if (pi != pj){
            //double r_2 = std::norm((x[pi] - x[pj], y[pi] - y[pj]));
            double r_2 = (pow(x[pi] - x[pj],2)+pow(y[pi] - y[pj],2));
            double ai = gravity_const * m[pj] / (r_2 + err);
            double ai_x = ai / (sqrt(r_2) + err) * (- x[pi] + x[pj]);
            double ai_y = ai / (sqrt(r_2) + err) * (- y[pi] + y[pj]);
            // std::cout << ((pow(ai_x,2) + pow(ai_y,2))) << std::endl;
            // std::cout << (pow(ai,2)) << std::endl;
            // std::cout << "ai_x " << ai_x << std::endl;
            // std::cout << "ai_y " << ai_y << std::endl;

            
            vx[pi] += ai_x * dt;
            vy[pi] += ai_y * dt;
            // std::cout << "vi_x " << vx[pi] << std::endl;
            // std::cout << "vi_y " << vy[pi] << std::endl;
            // std::cout << "-----------------------" << std::endl;
            if (sqrt(r_2) <= 2 * sqrt(radius2)){// collision
                double vec_r_xi = -x[pi] + x[pj];
                double vec_r_yi = -y[pi] + y[pj];

                double vec_v_xi = x[pi] + vx[pi];
                double vec_v_yi = y[pi] + vy[pi];
                double ci_2 = pow(vec_v_xi-x[pj],2)+pow(vec_v_yi-y[pj],2);
                double vi = sqrt(pow(vx[pi],2)+pow(vy[pi],2));
                double cos_ai =  (- ci_2 + r_2 + pow(vi,2)) / (2 * sqrt(r_2) * vi + err);
                double vi_changed = vi * cos_ai;
                //double vi_changed_ = vi * sqrt(1-pow(cos_ai,2));

                double vi_changed_x = - vi_changed * vec_r_xi / (vi + err);
                double vi_changed_y = - vi_changed * vec_r_yi / (vi + err);
                double vi_changed_x_ = vx[pi] + vi_changed_x;
                double vi_changed_y_ = vy[pi] + vi_changed_y;

                vx[pi] = vi_changed_x_ + vi_changed_x;
                vy[pi] = vi_changed_y_ + vi_changed_y;
            }

        }
    }

    // wall collision
    if (x[pi] >= 2000.0f + bound_x / 4 - 1){
        vx[pi] -= 2 * ori_vx;
    }
    if (y[pi] >= 2000.0f + bound_y / 4 - 1){
        vy[pi] -= 2 * ori_vy;
    }
    if (x[pi] <= 2000.0f - bound_x / 4 + 1){
        vx[pi] -= 2 * ori_vx;
    }
    if (y[pi] <= 2000.0f - bound_y / 4 + 1){
        vy[pi] -= 2 * ori_vy;
    }

}


void master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    //Logger l = Logger("openmp", n_body, bound_x, bound_y);
    double t_span = 0;
    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        
        //TODO: choose better threads configuration
        omp_set_num_threads(n_omp_threads);
        #pragma omp parallel for
        for (int i = 0; i < n_body; i++) {
            //std::cout << i << std::endl;
            update_velocity(m, x, y, vx, vy, i, n_body);
        }

        omp_set_num_threads(n_omp_threads);
        #pragma omp parallel for
        for (int i = 0; i < n_body; i++) {
            //std::cout << i << std::endl;
            update_position(x, y, vx, vy, i);
        }

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;
        t_span += time_span.count();
        //printf("Iteration %d, elapsed time: %.3f\n", i, time_span);

        //l.save_frame(x, y);

        #ifdef GUI
        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(2.0f);
        glBegin(GL_POINTS);
        double xi;
        double yi;
        for (int i = 0; i < n_body; i++){
            xi = x[i];
            yi = y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
        #else

        #endif
    }

    delete[] m;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;

    printf("%f ",t_span/n_iteration);

}


int main(int argc, char *argv[]){
    
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_omp_threads = atoi(argv[3]);

    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N Body Simulation Sequential Implementation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
	// if (n_omp_threads == 1)
	// 	printf("%d & ", n_body);

    master();

    // if (n_omp_threads == 40)
    //     printf("\\ \n");
    // else
    //     printf("& ");  

    printf("Student ID: 119010177\n"); // replace it with your student id
    printf("Name: Muhan Lin\n"); // replace it with your name
    printf("Assignment 2: N Body Simulation OpenMP Implementation\n");

    return 0;

}


