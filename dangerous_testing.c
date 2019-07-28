// principal code by Jasper Toscani Field
// with assistance from Tingjian Zhang

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#define KMER_SIZE 14
#define READ_LEN 300
#define READ_NUMBER 5001
#define LOCI_COUNT 10
#define LOCI_LEN 1000000
#define READ_KMER_NUM 30
//#define N_ROW 5001
//#define N_COL 2
//#define N_CHAR 300

#define send_data_tag 2001
#define first_send_data_tag 201
#define return_data_tag 2002

int main ( int argc, char *argv[] ) {

  FILE *fptr;
  FILE *msaptr;
  fptr = fopen(argv[1], "r"); // "r" for read
  msaptr = fopen(argv[2], "r");
  char read_count_buf[400];
  //char msa[LOCI_COUNT][1][LOCI_LEN];
  char read_buf[400];
  //char data2[READ_NUMBER][300];
  char msa_buf[LOCI_LEN];
  int i = 0;
  int j,k,x,n,z,loop;
  z = 0;
  int read_count = 0;
  int num_part = 0;
  int loci_count = 0;
  int total_read_count = 0;
  int total_seq_count = 0;
  int match_count = 0;
  const int STEPSIZE = 1000;
  int arrlen = STEPSIZE;
  const int MSA_STEPSIZE = 1;
  int msa_arrlen = MSA_STEPSIZE;
  const int MATCH_STEPSIZE = 1000;
  int match_arrlen = MATCH_STEPSIZE;
  const int MAX_ARRAY_SIZE = 1000000000;

  MPI_Status status;
  int my_id, root_process, ierr, num_rows, num_procs,
         an_id, num_rows_to_receive, avg_rows_per_process, 
         sender, num_rows_received, start_row, end_row, num_rows_to_send, avg_reads_per_process, avg_reads_to_receive, reads_to_receive;
  ierr = MPI_Init(&argc, &argv);
  root_process = 0;
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  char data_local[1000000];
  int data_pos_local[100000];


  if(my_id == root_process)
  {




// ARRAYS FOR HOLDING READ AND SEQUENCE STARTING POSITIONS IN BOTH OF THE FOLLOWING ARRAYS
//##############################################################################################
  int *data_pos = (int *)malloc(10000000 * sizeof(int));

  int *msa_pos = (int *)malloc(10000000 * sizeof(int));

//  char *data_local = (char *)malloc(1000000 * sizeof(char));

//  int *data_pos_local = (int *)malloc(100000 * sizeof(int));
//#############################################################################################

// CONSTRUCT A DYNAMIC MEMORY ARRAY TO HOLD THE READS WITHOUT REQUIRING KNOWING HOW MANY READS ARE IN THE FILE
//##############################################################################################  
  char *data = (char *)malloc(MAX_ARRAY_SIZE * sizeof(char));


  // READ READS FILE IN TO ARRAY
  while(fgets(read_buf, sizeof(read_buf), fptr))
  {
    
    
    read_count++;
    if(read_count == 2)
    {
      int start_pos = (strlen(data));
      if(start_pos == 0)
      {
	data_pos[total_read_count] = (start_pos);
      }
      else if(start_pos > 0)
      {
	data_pos[total_read_count] = (start_pos + 1);
      }

      //data_pos[total_read_count] = start_pos;
      //printf("%d\n", start_pos);
      read_buf[strlen(read_buf) - 1] = '\0';
      strcat(data, read_buf);
      
    }
    else if(read_count == 4)
    {
      read_count = 0;
      
      total_read_count++;
    }
  }

   /*
   printf("total read count %d\n", total_read_count);
   for(i = 0; i < strlen(data); i++)
   {
      printf("%c", data[i]);
   }
   */
  /*
  for(i = 0; i < total_read_count; i++)
  {
    printf("%d\n", data_pos[i]);
  }
  */
//##############################################################################################
  printf("\n\n\n DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\n\n\n");

// READ IN MSA FILE TO ARRAY
//#############################################################################################
  char *msa = (char *)malloc(MAX_ARRAY_SIZE * sizeof(char));
  //char **msa = (char **)malloc(MSA_STEPSIZE * sizeof(char *));

// READ IN MSA FILE TO ARRAY
  k = 0;
  while(fgets(msa_buf, sizeof(msa_buf), msaptr))
  {

    k++;
    if (k == 2)
    {
        msa_buf[strlen(msa_buf) - 1] = '\0';
	//printf("\n\n\n\n\n");
        int start_pos = (strlen(msa));
        if(start_pos == 0)
        {
          msa_pos[total_seq_count] = (start_pos);
        }
        else if(start_pos > 0)
        {
          msa_pos[total_seq_count] = (start_pos + 1);
        }
        int seq_len = strlen(msa_buf);
        //printf("seq leng = %d\n", seq_len);
        //char *seq = (char *)malloc((seq_len + 1) * sizeof(char));
        strcat(msa, msa_buf);

        //msa[total_seq_count] = seq;
        total_seq_count++;
        k = 0;

    }
   }
  
  //printf("total seq count %d\n", total_seq_count);
  /*
   for(i = 0; i < strlen(msa); i++)
   {
           printf("%c", msa[i]);
   }
  */
  /* 
  for(i = 0; i < total_seq_count; i++)
  {
    printf("%d\n", msa_pos[i]);
  }
  */
 
  avg_reads_per_process = total_read_count / num_procs;
  printf("average number of reads per processor = %d\n", avg_reads_per_process); 
  int sent_reads = 0;
  for(an_id = 1; an_id < num_procs; an_id++)
  {
    /*
    start_row = an_id*avg_rows_per_process + 1;
    end_row   = (an_id + 1)*avg_rows_per_process;

    if((total_read_count - end_row) < avg_rows_per_process)
      end_row = total_read_count - 1;

    num_rows_to_send = end_row - start_row + 1;
    */
    //sent_reads+=avg_reads_per_process;
    printf("%d\n", sent_reads); 
    int reads_to_send = READ_LEN * avg_reads_per_process;

    //SEND THE POSITION DATA ARRAY
    ierr = MPI_Send( &avg_reads_per_process, 1 , MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);

    ierr = MPI_Send( &data_pos[sent_reads], avg_reads_per_process, MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);


    //SEND THE STRING OF READS
    
    ierr = MPI_Send( &reads_to_send, 1 , MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);

    ierr = MPI_Send( &data[sent_reads], reads_to_send, MPI_CHAR,
           an_id, send_data_tag, MPI_COMM_WORLD);
    
    sent_reads+=avg_reads_per_process;
    printf("%d\n", sent_reads);

   }    
  }
//##########################################################################
//STUFF SENT TO OTHER PROCESSORS  
  else
  {
    //printf("PPPPPPPPPPPPPPPPPp%d\n", data_pos_local[1]);

    ierr = MPI_Recv( &avg_reads_to_receive, 1, MPI_INT, 
               root_process, send_data_tag, MPI_COMM_WORLD, &status);
    printf("wwwwwwwww%d\n", avg_reads_to_receive);      
    ierr = MPI_Recv( &data_pos_local, avg_reads_to_receive, MPI_INT, 
               root_process, send_data_tag, MPI_COMM_WORLD, &status);

    
    ierr = MPI_Recv( &reads_to_receive, 1, MPI_INT,
               root_process, send_data_tag, MPI_COMM_WORLD, &status);

    ierr = MPI_Recv( &data_local, reads_to_receive, MPI_CHAR,
               root_process, send_data_tag, MPI_COMM_WORLD, &status);
    
   printf("PPPPPPPPPPPPPPPPPp%c\n", data_local[1]);
   /*
   int reads_received = avg_reads_to_receive;
   for(i = 0; i < reads_received; i++)
  {
    printf("%d\n", data_pos_local[i]);
  }    
  */



  }
//END OF INDIVIDUAL PROCESSOR/CORE JOBS
//############################################################################



ierr = MPI_Finalize();

  return 0;
}




