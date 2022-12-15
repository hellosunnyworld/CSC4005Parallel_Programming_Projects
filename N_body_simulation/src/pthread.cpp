#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <pthread.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"

int n_thd; // number of threads

int n_body;
int n_iteration;


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

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n) {
    //TODO: calculate force and acceleration, update velocity

}


typedef struct {
    //TODO: specify your arguments for threads
    double* m;
    double* x;
    double* y;
    double* vx;
    double* vy;
    int n;
    int my_num;
    int tid;
    //TODO END
} Args;


void* worker(void* args) {
    //TODO: procedure in each threads

    Args* my_arg = (Args*) args;
    double* m = my_arg->m;
    double* x = my_arg->x;
    double* y = my_arg->y;
    double* vx = my_arg->vx;
    double* vy = my_arg->vy;
    int n = my_arg->n;
    int my_num = my_arg->my_num;
    int tid = my_arg->tid;
    double ori_vx;
    double ori_vy;
    
    for (int pi = tid * my_num; pi < tid * my_num + my_num; pi++){
        ori_vx = vx[pi];
        ori_vy = vy[pi];
        for (int pj = 0; pj < n; pj++){
            if (pi != pj){
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
                // std::cout << "vi_x " << vx[id] << std::endl;
                // std::cout << "vi_y " << vy[id] << std::endl;
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
    return NULL;
    // TODO END
}
void update_position(double x[],double y[],double vx[],double vy[], int n) {
    //TODO: update position 
    for (int i = 0; i < n; i++){
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
}

void master(){
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    //Logger l = Logger("pthread", n_body, bound_x, bound_y);
    int my_num = n_body / n_thd;
    double t_span = 0;
    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        //TODO: assign jobs
        pthread_t thds[n_thd]; // thread pool
        Args args[n_thd]; // arguments for all threads
        for (int thd = 0; thd < n_thd; thd++){
            args[thd].m = m;
            args[thd].x = x;
            args[thd].y = y;
            args[thd].vx = vx;
            args[thd].vy = vy;
            args[thd].n = n_body;
            args[thd].my_num = my_num;
            args[thd].tid = thd;
        }
        
        for (int thd = 0; thd < n_thd; thd++) pthread_create(&thds[thd], NULL, worker, &args[thd]);
        for (int thd = 0; thd < n_thd; thd++) pthread_join(thds[thd], NULL);
        int remainder = n_body - my_num * n_thd;
        if (remainder != 0){
            Args arg;
            arg.m = m;
            arg.x = x;
            arg.y = y;
            arg.vx = vx;
            arg.vy = vy;
            arg.n = n_body;
            arg.my_num = remainder;
            arg.tid = n_thd;
            worker(&arg);
        }
        update_position(x, y, vx, vy, n_body);
        //TODO End

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
        //pthread_exit(NULL);
    }
    // t_span /= n_iteration;
	// if (n_thd == 1)
	// 	printf("%d & ", n_body);
    // printf("%f ",t_span);
    // if (n_thd == 40)
    //     printf("\\ \n");
    // else
    //     printf("& ");
    delete[] m;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;


}


int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_thd = atoi(argv[3]);

    #ifdef GUI
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Pthread");
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    master();

	return 0;
}

