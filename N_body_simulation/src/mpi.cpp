#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"


int n_body;
int n_iteration;


int my_rank;
int world_size;


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


void update_position(double x[],double y[],double vx[],double vy[], int start, int end, double resultx[], double resulty[]) {
    //TODO: update position 
    for (int i = start; i < end; i++){
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

        resultx[i - start] = x[i];
        resulty[i - start] = y[i];
    }
}


void update_velocity(double m[], double x[],double y[],double vx[],double vy[], int my_num) {
    //TODO: calculate force and acceleration, update velocity
    double ori_vx;
    double ori_vy;
    int start, end;
    if (my_num == n_body / world_size){
        start = my_rank * my_num;
        end = my_rank * my_num + my_num;
    }
    else{
        start = world_size * my_num;
        end = n_body;
    }
    for (int pi = start; pi < end; pi++){
        ori_vx = vx[pi];
        ori_vy = vy[pi];
        for (int pj = 0; pj < n_body; pj++){
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
}


// void slave(int my_num){
//     // TODO: MPI routine
//     // TODO End

//     update_velocity(m, x, y, vx, vy, my_num);
// }



void master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    //Logger l = Logger("mpi", n_body, bound_x, bound_y);
    int my_num = n_body / world_size;
    MPI_Bcast(m, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vx, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vy, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double t_span = 0;
    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        MPI_Bcast(x, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(y, n_body, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // TODO: MPI routine

		if (my_rank == 0) {
			// you may change this part
			update_velocity(m,x,y,vx,vy,my_num);
		} 
		else {
			// you may change this part
			update_velocity(m,x,y,vx,vy,my_num);
		}     
        // TODO End
        MPI_Barrier(MPI_COMM_WORLD);

        double resultx[my_num];
        double resulty[my_num];
        if (my_rank != 0){
            update_position(x, y, vx, vy, my_rank * my_num, my_rank * my_num + my_num, resultx, resulty);
        }else{
			int remainder = n_body - my_num * world_size;
            double resultx0[my_num];
            double resulty0[my_num];
			if (remainder != 0){
				update_velocity(m,x,y,vx,vy,remainder);
                update_position(x, y, vx, vy,  world_size * my_num, n_body, resultx0, resulty0);
            }
            update_position(x, y, vx, vy, my_rank * my_num, my_rank * my_num + my_num, resultx, resulty);
        }
        
	    MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(resultx, my_num, MPI_DOUBLE, x, my_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(resulty, my_num, MPI_DOUBLE, y, my_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (my_rank == 0){
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
    }
    // if (my_rank == 0){
    //     if (world_size == 1)
    //         printf("%d & ", n_body);
    //     printf("%f ",t_span / n_iteration);
    //     if (world_size == 40)
    //         printf("\\ \n");
    //     else
    //         printf("& ");
    // }


    delete[] m;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;

}




int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	if (my_rank == 0) {
		#ifdef GUI
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(500, 500); 
		glutInitWindowPosition(0, 0);
		glutCreateWindow("N Body Simulation MPI Implementation");
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glMatrixMode(GL_PROJECTION);
		gluOrtho2D(0, bound_x, 0, bound_y);
		#endif
	} 
    master();

	if (my_rank == 0){
		printf("Student ID: 119010001\n"); // replace it with your student id
		printf("Name: Your Name\n"); // replace it with your name
		printf("Assignment 2: N Body Simulation MPI Implementation\n");
	}

	MPI_Finalize();

	return 0;
}

