#define _GNU_SOURCE
#include <sched.h>   
#include <unistd.h> 
#include <pthread.h>  
#include <math.h>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW  50000.0
#define ALIENS_HIGH 150000.0

signal* gloabl_sig;
double global_bandwidth;
int global_filter_order;
double* global_band_power;
int global_num_bands;
int global_num_threads;
int global_num_processors;

struct threadargs{
    long myid;
};

void usage() {
  printf("usage: p_band_scan text|bin|mmap signal_file Fs filter_order num_bands number-of-threads number-of-procs\n");
}


//Worker to make filter and run convolution
// Function run by each thread
void* worker(void* arg) {

    long myid     = (long)arg;
    signal* sig = gloabl_sig;
    double filter_coeffs[global_filter_order + 1];
    double bandwidth = global_bandwidth;
    int filter_order = global_filter_order;
    double* band_power = global_band_power;
    int num_bands = global_num_bands;
    int num_threads = global_num_threads;
    int  num_processors = global_num_processors;

    //divide between threads
    int blocksize = (global_num_bands + global_num_threads -1) / global_num_threads; // ceiling so last one might have fewer but more even load
    int mystart = myid * blocksize;
    int myend   = 0;

    if(mystart >= num_bands){ //we have more threads than we need so kill the thread
        pthread_exit(NULL);  
    }
    else if (myid == (num_threads - 1)) { // last thread
        myend = num_bands;
    } 
    else {
        myend = (myid + 1) * blocksize;
    }

    // put ourselves on the desired processor
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(myid % num_processors, &set);
    if (sched_setaffinity(0, sizeof(set), &set) < 0) { // do it
        perror("Can't setaffinity"); // hopefully doesn't fail
        exit(-1);
    }
    //calc for assigned bands
    for(int band = mystart; band<myend; band++){
        generate_band_pass(sig->Fs,
                            band * bandwidth + 0.0001, // keep within limits
                            (band + 1) * bandwidth - 0.0001,
                            filter_order,
                            filter_coeffs);
        hamming_window(filter_order,filter_coeffs);

        // Convolve
        convolve_and_compute_power(sig->num_samples,
                                    sig->data,
                                    filter_order,
                                    filter_coeffs,
                                    &(band_power[band]));
        }

    // Done.  The master thread will do the next thing
    pthread_exit(NULL);           // finish - no return value
}

void threadrunna(signal* sig, double bandwidth, int filter_order, double* band_power,int num_bands, int num_threads, int num_processors){
    pthread_t* tid = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);
    
    //assign gloabl vars so that threads can use em
    gloabl_sig = sig;
    global_filter_order = filter_order;
    global_band_power = band_power;
    global_num_bands = num_bands;
    global_num_threads = num_threads;
    global_num_processors = num_processors;
    global_bandwidth = bandwidth;
    
    if (!tid) {
        fprintf(stderr, "cannot allocate memory\n");
        exit(-1);
    }
    // launch threads
    for (long i = 0; i < num_threads; i++) {
        //call thread
        int returncode = pthread_create(&(tid[i]),  // thread id gets put here
                                        NULL, // use default attributes
                                        worker, // thread will begin in this function
                                        (void*)i // we'll give it i as the argument
                                        );
        if (returncode != 0) {
        perror("Failed to start thread");
        exit(-1);
        }
    }

    // now we will join all the threads
    for (int i = 0; i < num_threads; i++) {
        int returncode = pthread_join(tid[i], NULL);
        if (returncode != 0) {
        perror("join failed");
        exit(-1);
        }
    }
    free(tid);
}

double avg_power(double* data, int num) {

  double ss = 0;
  for (int i = 0; i < num; i++) {
    ss += data[i] * data[i];
  }

  return ss / num;
}

double max_of(double* data, int num) {

  double m = data[0];
  for (int i = 1; i < num; i++) {
    if (data[i] > m) {
      m = data[i];
    }
  }
  return m;
}

double avg_of(double* data, int num) {

  double s = 0;
  for (int i = 0; i < num; i++) {
    s += data[i];
  }
  return s / num;
}

void remove_dc(double* data, int num) {

  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (int i = 0; i < num; i++) {
    data[i] -= dc;
  }
}


int analyze_signal(signal* sig, int filter_order, int num_bands, double* lb, double* ub, int num_threads, int num_processors) {
  double Fc        = (sig->Fs) / 2;
  double bandwidth = Fc / num_bands;

  remove_dc(sig->data,sig->num_samples);

  double signal_power = avg_power(sig->data,sig->num_samples);

  printf("signal average power:     %lf\n", signal_power);

  resources rstart;
  get_resources(&rstart,THIS_PROCESS);
  double start = get_seconds();
  unsigned long long tstart = get_cycle_count();
  printf("Up to here took until %lld cycles, %f seconds", tstart, start);

//   double filter_coeffs[filter_order + 1];
  double band_power[num_bands];
//   for (int band = 0; band < num_bands; band++) {
//     // Make the filter
//     generate_band_pass(sig->Fs,
//                        band * bandwidth + 0.0001, // keep within limits
//                        (band + 1) * bandwidth - 0.0001,
//                        filter_order,
//                        filter_coeffs);
//     hamming_window(filter_order,filter_coeffs);

//     // Convolve
//     convolve_and_compute_power(sig->num_samples,
//                                sig->data,
//                                filter_order,
//                                filter_coeffs,
//                                &(band_power[band]));

//   }
  threadrunna(sig,bandwidth,filter_order,band_power,num_bands,num_threads,num_processors);
  
  unsigned long long tend = get_cycle_count();
  double end = get_seconds();
  printf("Up to here took until %lld cycles, %f seconds", tend, end);


  resources rend;
  get_resources(&rend,THIS_PROCESS);

  resources rdiff;
  get_resources_diff(&rstart, &rend, &rdiff);

  // Pretty print results
  double max_band_power = max_of(band_power,num_bands);
  double avg_band_power = avg_of(band_power,num_bands);
  int wow = 0;
  *lb = -1;
  *ub = -1;

  for (int band = 0; band < num_bands; band++) {
    double band_low  = band * bandwidth + 0.0001;
    double band_high = (band + 1) * bandwidth - 0.0001;

    printf("%5d %20lf to %20lf Hz: %20lf ",
           band, band_low, band_high, band_power[band]);

    for (int i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++) {
      printf("*");
    }

    if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
        (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {

      // band of interest
      if (band_power[band] > THRESHOLD * avg_band_power) {
        printf("(WOW)");
        wow = 1;
        if (*lb < 0) {
          *lb = band * bandwidth + 0.0001;
        }
        *ub = (band + 1) * bandwidth - 0.0001;
      } else {
        printf("(meh)");
      }
    } else {
      printf("(meh)");
    }

    printf("\n");
  }

  printf("Resource usages:\n\
User time        %lf seconds\n\
System time      %lf seconds\n\
Page faults      %ld\n\
Page swaps       %ld\n\
Blocks of I/O    %ld\n\
Signals caught   %ld\n\
Context switches %ld\n",
         rdiff.usertime,
         rdiff.systime,
         rdiff.pagefaults,
         rdiff.pageswaps,
         rdiff.ioblocks,
         rdiff.sigs,
         rdiff.contextswitches);

  printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\n"
         "Note that cycle count only makes sense if the thread stayed on one core\n",
         tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
  printf("Analysis took %lf seconds by basic timing\n", end - start);

  return wow;
}

int main(int argc, char* argv[]) {

  if (argc != 8) {
    usage();
    return -1;
  }

  char sig_type    = toupper(argv[1][0]);
  char* sig_file   = argv[2];
  double Fs        = atof(argv[3]);
  int filter_order = atoi(argv[4]);
  int num_bands    = atoi(argv[5]);
  int num_threads  = atoi(argv[6]);
  int num_procs    = atoi(argv[7]);

  assert(Fs > 0.0);
  assert(filter_order > 0 && !(filter_order & 0x1));
  assert(num_bands > 0);

  printf("type:     %s\n\
file:     %s\n\
Fs:       %lf Hz\n\
order:    %d\n\
bands:    %d\n",
         sig_type == 'T' ? "Text" : (sig_type == 'B' ? "Binary" : (sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE")),
         sig_file,
         Fs,
         filter_order,
         num_bands);

  printf("Load or map file\n");

  signal* sig;
  switch (sig_type) {
    case 'T':
      sig = load_text_format_signal(sig_file);
      break;

    case 'B':
      sig = load_binary_format_signal(sig_file);
      break;

    case 'M':
      sig = map_binary_format_signal(sig_file);
      break;

    default:
      printf("Unknown signal type\n");
      return -1;
  }

  if (!sig) {
    printf("Unable to load or map file\n");
    return -1;
  }

  sig->Fs = Fs;

  double start = 0;
  double end   = 0;
  if (analyze_signal(sig, filter_order, num_bands, &start, &end, num_threads, num_procs)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}

