#include <pthread.h>
#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cmath>
#include <mpi.h>

#define L 10 // weight parameter
#define NUM_LISTS 3 // number of lists
#define NUM_TASKS 100   // number of tasks per process
#define NUM_TASKS_TO_SHARE 10   // border of tasks-num above which we will share them as extra

// threads stuff
pthread_t threads[2];  // thread-descriptors: recv & exec
pthread_mutex_t mutex;  // mutex

int const STOP_RECV_MARKER = -1;   // stop-receive marker
int *tasks; // tasks array

int numProcesses;   // number of processes
int rank;   // rank of current process
int numRemainingTasks;  // number of remaining tasks
int numExecutedTasks;   // number of tasks executed during one iteration

double globalRes = 0;   // global result
double globalResSum;    // global result sum

// will we print only iteration summary or all info?
bool PRINT_SUMMARY_ONLY = true;

void initTasks(int *tasksArray, int numTasks, int iterCounter) {
    // initializing tasks weights
    for (int i = 0; i < numTasks; i++) {
        tasksArray[i] = abs(50 - i % 100) * abs(rank - (iterCounter % numProcesses)) * L;
    }
}

void executeTasks(const int *tasksArray) {
    for (int i = 0; i < numRemainingTasks; i++) {
        // blocking tasks queue
        pthread_mutex_lock(&mutex);

        // saving current task's weight in local variable
        int currTaskWeight = tasksArray[i];

        // unlocking tasks queue
        pthread_mutex_unlock(&mutex);

        for (int j = 0; j < currTaskWeight; j++) {
            globalRes += sin(j);    // increasing global sum
            usleep(tasksArray[i]);  // sleeping task-weight microseconds
        }

        numExecutedTasks++; // increasing number of executed tasks
    }

    numRemainingTasks = 0;  // no tasks remain (all executed)
}

void* executeStartRoutine(void* args) {
    // allocating memory for tasks array
    tasks = new int [NUM_TASKS];

    // times
    double startIterTime, iterTime = 0, minIterTime = 0, maxIterTime = 0;

    MPI_Status status;

    for (int i = 0; i < NUM_LISTS; i++) {
        // initializing tasks array
        initTasks(tasks, NUM_TASKS, i);

        numRemainingTasks = NUM_TASKS;
        numExecutedTasks = 0;
        int numExtraTasks = 0;

        startIterTime = MPI_Wtime();

        // executing our own tasks
        if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: EXEC :: EXECUTING OUR OWN TASKS..." << std::endl;
        executeTasks(tasks);
        if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: EXEC :: OUR OWN TASKS (" << numExecutedTasks << ") EXECUTED SUCCESSFULLY" << std::endl;

        for (int j = 0; j < numProcesses; j++) {
            // we finished our own tasks => sending requests for extra ones to all processes...
            if (j != rank) {
                // ...except of us
                if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: EXEC :: SENDING REQUEST FOR EXTRA TASKS TO PROC " << j << "..." << std::endl;
                MPI_Send(&rank, 1, MPI_INT, j, 0, MPI_COMM_WORLD);

                // receiving number of extra tasks
                if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: EXEC :: RECEIVING NUMBER OF EXTRA TASKS FROM PROC " << j << "..." << std::endl;
                MPI_Recv(&numExtraTasks, 1, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
                if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: EXEC :: NUMBER OF EXTRA TASKS FROM PROC " << j << " IS " << numExtraTasks << std::endl;

                if (numExtraTasks) {
                    // if there are extra tasks for us => receiving tasks themselves
                    if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: EXEC :: RECEIVING EXTRA TASKS FROM PROC " << j << "..." << std::endl;
                    MPI_Recv(tasks, numExtraTasks, MPI_INT, j, 1, MPI_COMM_WORLD, &status);
                    if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: EXEC :: RECEIVED EXTRA TASKS FROM PROC " << j << std::endl;

                    // executing extra tasks
                    if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: EXEC :: EXECUTING EXTRA TASKS..." << std::endl;
                    numRemainingTasks = numExtraTasks;  // extra tasks remain
                    executeTasks(tasks);    // executing extra tasks
                    if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: EXEC :: " << numExtraTasks << " EXTRA TASKS EXECUTED SUCCESSFULLY" << std::endl;
                }
            }
        }

        // elapsing time
        iterTime = MPI_Wtime() - startIterTime;

        // counting min & max iteration times
        MPI_Allreduce(&iterTime, &minIterTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&iterTime, &maxIterTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // all processes have to reach this point before going further
        MPI_Barrier(MPI_COMM_WORLD);

        // printing iteration summary for every process
        printf("%d :: EXECUTED TASKS:\t\t%d\n", rank, numExecutedTasks);
        printf("%d :: GLOBAL RESULT:\t\t%f\n", rank, globalRes);
        printf("%d :: ITERATION TIME:\t\t%f sec\n", rank, iterTime);
        printf("%d :: DIS-BALANCE:\t\t%f\n", rank, maxIterTime - minIterTime);
        printf("%d :: DIS-BALANCE (%%):\t\t%.3f%%\n", rank, (maxIterTime - minIterTime) / maxIterTime * 100);

        // summing global results
        MPI_Allreduce(&globalRes, &globalResSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // all processes have to reach this point before going further
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // sending STOP_RECV_MARKER == marker that we don't want to receive tasks anymore
    if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: EXEC :: SENDING STOP_RECV_MARKER..." << std::endl;
    MPI_Send(&STOP_RECV_MARKER, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);

    delete[] tasks;

    return NULL;
}

void* receiveStartRoutine(void* args) {
    int numShareTasks, requestRank;

    MPI_Status status;

    while (true) {
        // receiving request
        if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: RECV :: RECEIVING REQUEST-RANK..." << std::endl;
        MPI_Recv(&requestRank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: RECV :: RECEIVED REQUEST-RANK IS " << requestRank << std::endl;

        if (requestRank == STOP_RECV_MARKER) {
            // received STOP_RECV_MARKER => finish
            if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: RECV :: RECEIVED STOP_RECV_MARKER -> EXIT" << std::endl;
            pthread_exit(NULL);
        }

        // blocking tasks queue
        pthread_mutex_lock(&mutex);

        if (numRemainingTasks > NUM_TASKS_TO_SHARE) {
            // we got a lot of tasks => let's share half of them as extra for asking proc
            numShareTasks = numRemainingTasks / 2;
            numRemainingTasks -= numShareTasks;

            // sending number of extra tasks & tasks themselves
            if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: RECV :: SENDING NUMBER OF EXTRA TASKS (" << numShareTasks << ")..." << std::endl;
            MPI_Send(&numShareTasks, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
            if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: RECV :: SENDING EXTRA TASKS THEMSELVES..." << std::endl;
            MPI_Send(&tasks[numShareTasks - 1], numShareTasks, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
        } else {
            // we don't want to share tasks

            numShareTasks = 0;

            if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: RECV :: SENDING NUMBER OF EXTRA TASKS (" << numShareTasks << ")..." << std::endl;
            MPI_Send(&numShareTasks, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
        }

        // unlocking tasks queue
        pthread_mutex_unlock(&mutex);
    }
}

int createThreads() {
    pthread_attr_t attrs;   // thread attributes

    // initializing thread attributes
    if (pthread_attr_init(&attrs)){
        perror("ERROR :: COULD NOT INITIALIZE ATTRIBUTES");
        return 1;
    }

    // setting 'joinable' for threads
    if (pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE)){
        perror("ERROR :: COULD NOT SET DETACH STATE");
        return 1;
    }

    // creating receive-thread
    if (pthread_create(&threads[0], &attrs, receiveStartRoutine, NULL)){
        perror("ERROR :: COULD NOT CREATE RECEIVE-THREAD");
        return 1;
    }

    // creating execute-thread
    if (pthread_create(&threads[1], &attrs, executeStartRoutine, NULL)){
        perror("ERROR :: COULD NOT CREATE EXECUTE-THREAD");
        return 1;
    }

    // setting attribute resources free
    pthread_attr_destroy(&attrs);

    // current thread doesn't stop until threads[0] finishes
    if (pthread_join(threads[0], NULL)) {
        perror("ERROR :: COULD NOT JOIN RECEIVE-THREAD");
        return 1;
    }

    // current thread stops until threads[1] finishes
    if (pthread_join(threads[1], NULL)){
        perror("ERROR :: COULD NOT JOIN EXECUTE-THREAD");
        return 1;
    }

    // nobody calls pthread_join for joinable thread => after finish it doesn't free its resources => memory leak

    return 0;
}

void tasksDistribution() {
    // getting number of processes and rank of current process
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: START" << std::endl;

    // initializing mutex before using
    // NUL = default attributes
    pthread_mutex_init(&mutex, NULL);

    // creating & starting threads
    createThreads();

    // destroying mutex after using
    pthread_mutex_destroy(&mutex);
}

int main(int argc, char **argv) {
    // initializing PRINT_SUMMARY_ONLY parameter with true for 1 or false for 0
    PRINT_SUMMARY_ONLY = (atoi(argv[1]) != 0);

    int provided;   // actually provided level of "thread support"

    // initializing multi-thread MPI
    // classic MPI_Init == MPI_Init_thread with flag MPI_THREAD_SINGLE
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE) {
        printf("ERROR :: COULD NOT GET NEEDED LEVEL\n");
        MPI_Finalize();
        return 1;
    }

    double startTime = MPI_Wtime();

    tasksDistribution();

    double elapsedTime = MPI_Wtime() - startTime;

    if (rank == 0) {
        std::cout << "ELAPSED TIME: " << elapsedTime << " sec" << std::endl;
        std::cout << "GLOBAL RESULT SUM: " << globalResSum << std::endl;
    }

    if (!PRINT_SUMMARY_ONLY) std::cout << rank << " :: FINISH" << std::endl;

    MPI_Finalize();

    return 0;
}