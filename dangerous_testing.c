// principal code by Jasper Toscani Field
// with assistance from Tingjian Zhang

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#define KMER_SIZE 14
#define READ_NUMBER 5001
#define LOCI_COUNT 10
#define LOCI_LEN 100
#define READ_KMER_NUM 30
//#define N_ROW 5001
//#define N_COL 2
//#define N_CHAR 300
#define send_data_tag 2001
#define return_data_tag 2002

int main ( int argc, char *argv[] ) {

  FILE *fptr;
  FILE *msaptr;
  fptr = fopen(argv[1], "r"); // "r" for read
  msaptr = fopen(argv[2], "r");
  char read_count_buf[400];
  //char msa[LOCI_COUNT][1][LOCI_LEN];
  char read_buf[400];
  char data2[READ_NUMBER][300];
  char msa_buf[LOCI_LEN];
  int i = 0;
  int j,k,x,n,z,loop;
  int read_count = 0;
  int num_part = 0;
  int loci_count = 0;
  int total_read_count = 0;
  int total_seq_count = 0;
  int match_count = 0;
  const int STEPSIZE = 100;
  int arrlen = STEPSIZE;
  const int MSA_STEPSIZE = 1;
  int msa_arrlen = MSA_STEPSIZE;
  const int MATCH_STEPSIZE = 1000;
  int match_arrlen = MATCH_STEPSIZE;
  MPI_Status status;
  int my_id, root_process, ierr, num_rows, num_procs,
         an_id, num_rows_to_receive, avg_rows_per_process, 
         sender, num_rows_received, start_row, end_row, num_rows_to_send;
  ierr = MPI_Init(&argc, &argv);
  root_process = 0;
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  if(my_id == root_process)
  {





// CONSTRUCT A DYNAMIC MEMORY ARRAY TO HOLD THE READS WITHOUT REQUIRING KNOWING HOW MANY READS ARE IN THE FILE
//##############################################################################################  
  char **data = (char **)malloc(STEPSIZE * sizeof(char *));


  // READ READS FILE IN TO ARRAY
  while(fgets(read_buf, sizeof(read_buf), fptr))
  {

    if ( z == arrlen)
    {
    arrlen += STEPSIZE;
    char ** newlines = realloc(data, arrlen * sizeof(char * ));
    if (!newlines)
    {
      fprintf(stderr, "can't realloc\n");
      exit(1);
    }
    data = newlines;
    }
    read_count++;
    if(read_count == 2)
    {
      //printf("%s\n", read_buf);
      read_buf[strlen(read_buf) - 1] = '\0';
      int read_len = strlen(read_buf);
      char *read_seq = (char *)malloc((read_len + 1) * sizeof(char));
      strcpy(read_seq, read_buf);

      data[z] = read_seq;

    }
    else if(read_count == 4)
    {
      read_count = 0;
      z++;
      total_read_count++;
    }
  }
   printf("total read count %d\n", total_read_count);
//##############################################################################################


// READ IN MSA FILE TO ARRAY
//#############################################################################################
  /*
   k = 0;
  while(fscanf(msaptr,"%s",msa[k][0]) != EOF)
  {
    fscanf(msaptr,"%s",msa[k][0]);
    k++;
  }
  */

  char **msa = (char **)malloc(MSA_STEPSIZE * sizeof(char *));

// READ IN MSA FILE TO ARRAY
  k = 0;
  while(fgets(msa_buf, sizeof(msa_buf), msaptr))
  {

    if ( total_seq_count == msa_arrlen)
    {
        msa_arrlen += MSA_STEPSIZE;
        char ** msa_newlines = realloc(msa, msa_arrlen * sizeof(char * ));
        if (!msa_newlines)
        {
                fprintf(stderr, "can't realloc\n");
                exit(1);
        }
        msa = msa_newlines;
        //printf("d\n");
    }
    k++;
    if (k == 2)
    {
        //printf("%s\n", msa_buf);
        //printf("k = %d\n", k);
        //printf("seq number = %d\n", total_seq_count);
        //printf("%s\n", msa_buf);
        msa_buf[strlen(msa_buf) - 1] = '\0';
        int seq_len = strlen(msa_buf);
        //printf("seq leng = %d\n", seq_len);
        char *seq = (char *)malloc((seq_len + 1) * sizeof(char));
        strcpy(seq, msa_buf);

        msa[total_seq_count] = seq;
        total_seq_count++;
        k = 0;

    }
   }


  printf("%s\n", data[0]);

  avg_rows_per_process = total_read_count / num_procs;
  printf("average number of reads per processor = %d\n", avg_rows_per_process);

  for(an_id = 1; an_id < num_procs; an_id++)
  {
    start_row = an_id*avg_rows_per_process + 1;
    end_row   = (an_id + 1)*avg_rows_per_process;

    if((total_read_count - end_row) < avg_rows_per_process)
      end_row = total_read_count - 1;

    num_rows_to_send = end_row - start_row + 1;

    ierr = MPI_Send( &num_rows_to_send, 1 , MPI_CHAR,
           an_id, send_data_tag, MPI_COMM_WORLD);

    ierr = MPI_Send( &data[start_row][0], num_rows_to_send, MPI_CHAR,
           an_id, send_data_tag, MPI_COMM_WORLD);
   }
  }
  else
  {
    ierr = MPI_Recv( &num_rows_to_receive, 1, MPI_INT, 
               root_process, send_data_tag, MPI_COMM_WORLD, &status);
          
    ierr = MPI_Recv( &data2, num_rows_to_receive, MPI_INT, 
               root_process, send_data_tag, MPI_COMM_WORLD, &status);

    num_rows_received = num_rows_to_receive;
    for(i = 0;i < num_rows_received; i++);
    {
      for(j = 0; j < 300; j++)
      {
	      printf("%s",data2[i]);
      }
      printf("\n");
    }





    printf("PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n");
    printf("Hello world! I'm process %i out of %i processes\n",
         my_id, num_procs);
  }
  
  printf("Finish Read\n");
  ierr = MPI_Finalize();

  return 0;
}
