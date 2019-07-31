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
#define send_data_tag 2001
#define first_send_data_tag 201
#define return_data_tag 2002

int main ( int argc, char *argv[] ) {

  FILE *fptr;
  FILE *msaptr;
  fptr = fopen(argv[1], "r"); // "r" for read
  msaptr = fopen(argv[2], "r");
  char read_count_buf[400];
  char read_buf[400];
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
  int avg_read_size;
  int reads_nuc_sum;
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
         sender, num_rows_received, start_row, end_row, num_rows_to_send, avg_reads_per_process, avg_reads_to_receive, reads_to_receive,
	 start, stop, msa_nuc_sum, avg_seq_len, avg_seq_len_local, total_seq_count_local;
  ierr = MPI_Init(&argc, &argv);
  root_process = 0;
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  char data_local[1000000];
  int data_pos_local[1000000];
  //int data_pos_local_buffer[1000000];
  char msa_data_local[1000000];
  char read_buffer[400];
  int match[1000000];
  if(my_id == root_process)
  {




// ARRAYS FOR HOLDING READ AND SEQUENCE STARTING POSITIONS IN BOTH OF THE FOLLOWING ARRAYS
//##############################################################################################
  int *data_pos = (int *)malloc(10000000 * sizeof(int));

  //int *msa_pos = (int *)malloc(10000000 * sizeof(int));

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
	reads_nuc_sum+=strlen(read_buf);
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
  avg_read_size = reads_nuc_sum / total_read_count;
  printf("all nucs in read seq = %d\n", reads_nuc_sum);
  /* 
  for(i = 0; i < strlen(data); i++)
  {
    printf("%c", data[i]);
  }
  */
  printf("\n");
  
//##############################################################################################

// READ IN MSA FILE TO ARRAY
//#############################################################################################
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
	msa_nuc_sum+=seq_len;
        //printf("seq leng = %d\n", seq_len);
        char *seq = (char *)malloc((seq_len + 1) * sizeof(char));
        strcpy(seq, msa_buf);

        msa[total_seq_count] = seq;
        total_seq_count++;
        k = 0;

    }
  }
  printf("FINISHED READING IN DATA FILES, BEGINNING SENDING DATA TO PROCESSORS\n"); 
  avg_seq_len = msa_nuc_sum / total_seq_count;
  //printf("%d\n", avg_seq_len);

 
  avg_reads_per_process = total_read_count / num_procs;
  printf("average number of reads per processor = %d\n", avg_reads_per_process); 
  int sent_reads = 0;
  int nucs_being_sent = 0;
  int starting_nuc = 0;
  printf("READ POSITIONS IN DATA_POS = \n");
  /*
  for( i = 0; i < total_read_count; i++)
  {
    printf("data_pos %d = %d\n", i, data_pos[i]);
  }
  */

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
    printf("the numbers of positional data sent = %d\n", sent_reads);
    printf(" the starting nucleotide being sent = %d\n", starting_nuc);
    //printf(" the number of nucleotides being sent = %d\n", reads_to_send); 
    int reads_to_send = avg_read_size * (avg_reads_per_process);
    printf(" the number of nucleotides being sent = %d\n", reads_to_send);
    printf("SENDING POSITIONAL DATA FOR START AND STOP OF READS\n");
    //SEND THE POSITION DATA ARRAY
    ierr = MPI_Send( &avg_reads_per_process, 1 , MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);

    ierr = MPI_Send( &data_pos[sent_reads], avg_reads_per_process, MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);
    printf("FINISHED SENDING POSITIONAL READ DATA\n");

    //SEND THE STRING OF READS
    printf("SENDING READ DATA\n");
    printf("SENDING THE NUMBER OF THE NUMBER OF READS TO SEND\n"); 
    ierr = MPI_Send( &reads_to_send, 1 , MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);
    printf("SENDING CONCATENATED READ SEQUENCE\n");
    //printf("reads to send = %d\n", reads_to_send);
    //printf("reads_to_send LS = %ls\n", &reads_to_send);
    //printf("sent reads = %d\n", sent_reads);
    //printf("length of read sequence = %d\n", strlen(data));
    //printf("length of read sequence LD = %ld\n", strlen(data));
    ierr = MPI_Send( &data[starting_nuc], reads_to_send, MPI_CHAR,
           an_id, send_data_tag, MPI_COMM_WORLD);
    
    sent_reads+=avg_reads_per_process;
    starting_nuc+=reads_to_send;
    printf("FINISHED SENDING READ DATA AND READ POSITIONAL INFO\n");


    //SEND MSA DATA
    printf("SENDING MSA DATA\n");
    ierr = MPI_Send( &avg_seq_len, 1 , MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);

    ierr = MPI_Send( &total_seq_count, 1 , MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);
    for(n = 0; n < total_seq_count; n++)
    {
      ierr = MPI_Send( &msa[n][0], avg_seq_len, MPI_CHAR,
           an_id, send_data_tag, MPI_COMM_WORLD);
    }

   printf("SENT ALL DATA TO A SINGLE PROCESSOR\n");
   }
  printf("FINISHED SENDING ALL INFO FROM PROCESS 1\n");  
  }
//##########################################################################
//STUFF SENT TO OTHER PROCESSORS  
  else
  {
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    //printf("PPPPPPPPPPPPPPPPPp%d   and ID = %i\n", data_pos_local[416611], my_id);
    //RECEIVE DATA FOR READS and READ START AND STOPPING LOCATIONS
    ierr = MPI_Recv( &avg_reads_to_receive, 1, MPI_INT, 
               root_process, send_data_tag, MPI_COMM_WORLD, &status);
    //printf("wwwwwwwww%d\n", avg_reads_to_receive);      
    ierr = MPI_Recv( &data_pos_local, avg_reads_to_receive, MPI_INT, 
               root_process, send_data_tag, MPI_COMM_WORLD, &status);
    
    //int data_pos_local_buffer[1000000];
    //data_pos_local_buffer[1] = 0;
    /* 
    for( i = 1; i < avg_reads_to_receive; i++)
    {
      //data_pos_local[i] = (data_pos_local_buffer[i + 1] - data_pos_local_buffer[i]) + data_pos_local_buffer[i - 1];
      printf("data pos local position %d =   data pos local value = %d\n", i, data_pos_local[i]);
    }
    */
    ierr = MPI_Recv( &reads_to_receive, 1, MPI_INT,
               root_process, send_data_tag, MPI_COMM_WORLD, &status);

    ierr = MPI_Recv( &data_local, reads_to_receive, MPI_CHAR,
    	    root_process, send_data_tag, MPI_COMM_WORLD, &status);

    /*
    int pos = 0;
    int counted_nucs = 0;
      for( z = 0; z < avg_reads_to_receive; z++)
      {
        loop = 0;
        if(data_pos_local[z + 1] != 0)
        {
          pos = (data_pos_local[z + 1] - data_pos_local[z]);
          printf("POS = %d, ID = %d\n", pos, my_id);
        }
        //printf("%c", data_local[i]);
        //printf("start = %d stop = %d\n",start , stop);
        //printf("%d\n", i); 
        for( j = 0; j < pos; j++)
        {
          //printf("%c", data_local[j]);
          read_buffer[j] = data_local[counted_nucs + j];
	  printf("%c", read_buffer[j]);
          //counted_nucs+=pos;
          counted_nucs++;
        }
	printf("\n");
      }
    */
    /*
    for(i = 0; i < avg_reads_to_receive; i++)
    {
      printf("ORIGINAL POSITION = %d, ID = %d\n", data_pos_local[i], my_id);
      if(data_pos_local[i + 1] != 0)
      {
	pos = (data_pos_local[i + 1] - data_pos_local[i]);
	printf("POS = %d, ID = %d\n", pos, my_id);
      }
        
      
    }
    */
    //printf("\n");  
    
    // RECEIVE DATA FOR MSA
    ierr = MPI_Recv( &avg_seq_len_local, 1, MPI_INT,
               root_process, send_data_tag, MPI_COMM_WORLD, &status);

    ierr = MPI_Recv( &total_seq_count_local, 1, MPI_INT,
               root_process, send_data_tag, MPI_COMM_WORLD, &status);
    
    
//###################################################################################
//RECEIVE MSA INDIVIDUAL SEQUENCES AND BEGIN LOOPING READS OVER THEM    
    for(n = 0; n < total_seq_count_local; n++)
    {
	    
      ierr = MPI_Recv( &msa_data_local, avg_seq_len_local, MPI_CHAR,
   	      root_process, send_data_tag, MPI_COMM_WORLD, &status);
       
      char *str = msa_data_local;
      //printf("%s\n", str);
      //char *read;
      //char read_buffer[400];
      char *kmer;
      char *nuc;
      int seq_size = strlen(str);
      //int kmer_size;
      int safe_hasher = 0;

      /* LOOP THROUGH READS AND BEGIN CUTTING INTO KMERS */
      printf("BEGIN LOOPING THROUGH READS\n");
      //int start_num = 0;
      //int stop_num = 1;
      //start = data_pos_local[start_num];
      //stop = data_pos_local[stop_num];
      //printf("START = %d STOP = %d", start, stop);
      int pos = 0;
      int counted_nucs = 0;
      for( z = 0; z < avg_reads_to_receive - 2; z++)
      {
	loop = 0;
        if(data_pos_local[z + 1] != 0)
        {
          pos = (data_pos_local[z + 1] - data_pos_local[z]);
          //printf("POS = %d, ID = %d\n", pos, my_id);
        }
        //printf("%c", data_local[i]);
        //printf("start = %d stop = %d\n",start , stop);
        //printf("%d\n", i); 
        for( j = 0; j < pos; j++)
        {
          //printf("%c", data_local[j]);
          read_buffer[j] = data_local[counted_nucs + j];
          //printf("%c", read_buffer[j]);
          //counted_nucs+=pos;
          counted_nucs++;
        }
        //printf("\n");
	//printf("\n##################################### READ BUFFER OUTPUT TO FOLLOW");
        for( loop = 0; loop < strlen(read_buffer); loop++)
	{
	  //printf("%c", read_buffer[loop]);
	}
        //printf("\n");
        int len = strlen(read_buffer);
        char *read = read_buffer;
        int read_size = strlen(read);
        int read_km_number = read_size / KMER_SIZE;
        int read_hash_array[1][read_km_number + 2];
        int trimmed_read_size = read_size - read_size % KMER_SIZE;
        int window = 0;
        for(j = 0; j < read_km_number; j++)
        {
          int read_hash_number = 0;
          for (i = window; i < KMER_SIZE * (j + 1); i++)
          {
            window++;
            int read_hash = 0;
            int ps = read[i];
            //printf("%c",ps); 
        
            /* HASHES GENERATED FOR THE KMERS FROM THE READ */
	    switch(ps)
	    {
	      case 'A':
	      case 'a':
	        read_hash = 0;
	        break;
	      case 'C':
              case 'c':
                read_hash = 1;
                break;
	      case 'G':
              case 'g':
                read_hash = 2;
                break;
	      case 'T':
              case 't':
                read_hash = 3;
                break;
	      default :
	        read_hash = 4;
	    }
            read_hash_number = read_hash_number * 4 + read_hash;
	    read_hash = 0;
          }
          read_hash_array[0][j] = read_hash_number;
          read_hash_number = 0;
        }


        //printf("***SEQUENCE HASHING SLIDING WINDOW***\n");
        int msa_hash_number = 0;
	printf("msa_hash_number = %d\n", msa_hash_number);
        for(i = 0 ; i < seq_size; i++)
        {
          printf("%c\n", str[i]);
          int hash = 0;
	  switch(str[i])
          {
            case 'A':
            case 'a':
              hash = 0;
              break;
            case 'C':
            case 'c':
              hash = 1;
              break;
            case 'G':
            case 'g':
              hash = 2;
              break;
            case 'T':
            case 't':
              hash = 3;
              break;
	    default :
	      hash = 4;
          }
	  int sub_hash = 0;
	  switch(str[i - KMER_SIZE])
          {
            case 'A':
            case 'a':
              sub_hash = 0;
              break;
            case 'C':
            case 'c':
              sub_hash = 1;
              break;
            case 'G':
            case 'g':
              sub_hash = 2;
              break;
            case 'T':
            case 't':
              sub_hash = 3;
              break;
	    default :
	      sub_hash = 4;
          }
          msa_hash_number = msa_hash_number * 4 + hash - (sub_hash << 28);
          for(j = 0; j < read_km_number; j++)
          {
            printf("read kmer = %d   msa kmer = %d\n", read_hash_array[0][j], msa_hash_number);
           if(read_hash_array[0][j] == msa_hash_number)
           {
//############################################################
              printf("FOUND HASH MATCH seq position:%d in loci: %d read kmer:%d of read %d\n", i, n, j, z);
              //printf("%s\n", read);
	      //match_count++;
	      /*
	      if(match_count == match_arrlen)
	      {
	        match_arrlen += MATCH_STEPSIZE;
                int ** match_newlines = realloc(match, match_arrlen * sizeof(int * ));
                if (!match_newlines)
                {
                  fprintf(stderr, "can't realloc\n");
                  exit(1);
                }
                match = match_newlines;
	      }
	      */
	      //match[match_count] = (int *)malloc(4 * sizeof(int));
	      //match[match_count] = x;
	      //match[match_count + 1] = j;
	      //match[match_count + 2] = n;
	      //match[match_count + 3] = i;
	      //match_count+=5;
	    
//###################################################################
//FINISH LOOPING OVER INDIVIDUAL NUCS OF MSA AND CONSTRUCTING HASH'S
	  }
        }
	printf("msa hash = %d\n", hash);
	hash = 0;
	sub_hash = 0;
//####################################################################

//####################################################################
      }
      msa_hash_number = 0;
//FINISH LOOPING OVER MSA AND DOING COMPARISON
//#####################################################################

//#####################################################################
//FINISH LOOPING OVER CONTIGUOUS READ STRING AND MOVE TO NEXT READ
    //printf("\n");
    //start_num++;
    //stop_num++;
    //start = data_pos_local[start_num];
    //stop = data_pos_local[stop_num];
    }
//#####################################################################

//#####################################################################
//FINISH LOOPING OVER ALL MSA SEQUENCES COMING FOR ROOT PROCESSOR
   }


//#####################################################################



 
  //free(data);
  /* 
  for(i = 0; i < match_count; i++)
  {
	  for(j = 0; j < 4; j++)
	  {
		  printf("  %d   ", match[i][j]);
		  //free(&match[i][j]);
	  }
	  printf("\n");
          free(match[i]);
  }
  free(match);
  */

//####################################################################
//FINISH ALL PROCESSES BY SUBSERVIENT PROCESSORS
    }
//####################################################################

   /*
   start_num = 0;
   stop_num = 1;
   start = data_pos_local[start_num];
   stop = data_pos_local[stop_num];
   //printf("START = %d STOP = %d", start, stop);
   for( i = 0; i < avg_reads_to_receive; i++)
   {
      loop = 0;
      //printf("%c", data_local[i]);
      //printf("start = %d stop = %d\n",start , stop);
     
     for( j = start; j < stop; j++)
     {
       //printf("%c", data_local[j]);
       read_buffer[loop] = data_local[j];
       loop++;
     }
     */






    /*

     //printf("\n");
     start_num++;
     stop_num++;
     start = data_pos_local[start_num];
     stop = data_pos_local[stop_num]; 
   }
   printf("\n%d\n", avg_reads_to_receive);
   */
//
//END OF INDIVIDUAL PROCESSOR/CORE JOBS
//############################################################################



ierr = MPI_Finalize();

  return 0;
}




