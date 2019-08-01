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
  int j,k,x,n,z,loop,p,q,r;
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
      int start_pos = (strlen(read_buf));
      //int start_pos = (strlen(data));
      
      data_pos[total_read_count] = start_pos;
      
      /*
      if(start_pos == 0)
      {
	data_pos[total_read_count] = (start_pos);
      }
      else if(start_pos > 0)
      {
	data_pos[total_read_count] = (start_pos + 1);
	reads_nuc_sum+=strlen(read_buf);
      }
      */
      reads_nuc_sum+=strlen(read_buf);
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
  printf("total read count for run is = %d\n", total_read_count);
  printf("all nucs in read seq = %d\n", reads_nuc_sum);
  /*
  for(i = 0; i < total_read_count; i++)
  {
	  printf("%d\n", data_pos[i]);
  }
  */
  //printf("\n");
  
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
    }
    k++;
    if (k == 2)
    {
        msa_buf[strlen(msa_buf) - 1] = '\0';
        int seq_len = strlen(msa_buf);
	msa_nuc_sum+=seq_len;
        char *seq = (char *)malloc((seq_len + 1) * sizeof(char));
        strcpy(seq, msa_buf);

        msa[total_seq_count] = seq;
        total_seq_count++;
        k = 0;

    }
  }
  printf("FINISHED READING IN DATA FILES, BEGINNING SENDING DATA TO PROCESSORS\n"); 
  avg_seq_len = msa_nuc_sum / total_seq_count;
  int track_leftover_nucs = 0; 
  //printf("%d\n", avg_seq_len);
  printf("BEGINNING CALCULATIONS FOR SENDING DATA\n");
  //int reads_to_send = avg_read_size * (avg_reads_per_process);
  avg_reads_per_process = total_read_count / num_procs;
  printf("Average number of reads per core = %d\n", avg_reads_per_process);
  int reads_to_send = avg_read_size * (avg_reads_per_process);
  printf("the number of nucleotides sent to each core will be = %d\n", reads_to_send);
  int leftover_reads = total_read_count % num_procs;
  //printf("the number of leftover reads = %d\n", leftover_reads);
  int leftover_reads_len = avg_read_size * leftover_reads;
  //printf("the length of the leftover reads is = %d nucleotides\n", leftover_reads_len);
  /*
  int sent_reads = 0;
  int nucs_starting_pos = 0;
  int reads_starting_pos = 0;
  int read_sending_counter = 0;
  int reads_this_stage = 0;
  printf("CALCULATIONS COMPLETE, SENDING DATA \n");
  int send_stage = 0;
  int an_id = 1;
  int read_count = avg_reads_per_process;
  */
  int numbers_of_reads[num_procs];
  int procs_read_numbers[num_procs + 1];
  procs_read_numbers[0] = 0;
  int procs_reads_starting_pos[num_procs + 1];
  procs_reads_starting_pos[0] = 0;
  int procs_nucleotide_number[num_procs + 1];
  //procs_nucleotide_number[0] = 0;
  int procs_read_start_position[num_procs + 1];
  procs_read_start_position[0] = 0;
  int procs_nucleotide_start_position[num_procs + 1];
  procs_nucleotide_start_position[0] = 0;
  int read_count = 0;
  int read_count_buffer = 0;
  int proc_count = 1;
  int read_len;
  int total_nucs = 0;
  int read_len_counter = 0;
  int reads_per_proc;
  my_id = 1;
  for(i = 0; i < total_read_count; i++)
  {
     read_len = data_pos[i];
     read_count += read_len;
     read_len_counter++;
     total_nucs += read_len;
     if(read_len_counter == avg_reads_per_process || i == (total_read_count - 1))
     {
       printf("nucleotides for each  proccessor = %d\n", read_count);
       //reads_per_proc = 
       procs_nucleotide_number[my_id] = read_count;
       procs_read_numbers[my_id] = i;
       procs_nucleotide_start_position[my_id] = total_nucs;
       printf("nucleotide starting position = %d\n", total_nucs);
       printf("i = %d\n", i);
       read_count = 0;
       read_len_counter = 0;
       my_id++;
     }
  }
   
  for(an_id = 1; an_id < num_procs; an_id++)
  {
        
    //nucs_starting_pos = reads_this_stage;
    //reads_starting_pos = read_count;
    //send info message
    ierr = MPI_Send( &procs_read_numbers[an_id], 1 , MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);

    ierr = MPI_Send( &data_pos[procs_read_numbers[an_id - 1]], procs_read_numbers[an_id], MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);
    printf("FINISHED SENDING POSITIONAL READ DATA\n");

         //SEND THE STRING OF READS
    ierr = MPI_Send( &procs_nucleotide_number[an_id], 1 , MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);

    ierr = MPI_Send( &data[procs_nucleotide_start_position[an_id - 1]], procs_nucleotide_number[an_id], MPI_CHAR,
           an_id, send_data_tag, MPI_COMM_WORLD); 
  
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






   }
  
  
  
  
  
  
  
  
  
  
  
  /*
  for(j = 0; j < total_read_count; j++)
  {
    sent_reads = data_pos[j];
    reads_this_stage += sent_reads;
    
    if((j == total_read_count - 1) || j >= read_count)
    {  

       if(an_id == num_procs - 1)
       {
         int last_reads = total_read_count - read_count;
         int last_nucs = reads_nuc_sum - nucs_starting_pos; 
         printf("leftover_reads to send = %d\n", last_reads);
         printf("leftover nucs = %d\n", last_nucs);
         
         //nucs_starting_pos = reads_this_stage;
         //reads_starting_pos = read_count;
         printf("FINAL reads this stage = %d     current step = %d    for id = %d\n", last_nucs, last_reads, an_id);
         printf("FINAL nucleotide starting position = %d     reads starting position = %d\n", nucs_starting_pos, reads_starting_pos);
         //send info message
         ierr = MPI_Send( &last_reads, 1 , MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);

         ierr = MPI_Send( &data_pos[reads_starting_pos], last_reads, MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);
         printf("FINISHED SENDING POSITIONAL READ DATA\n");

         //SEND THE STRING OF READS
         ierr = MPI_Send( &last_nucs, 1 , MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);

         ierr = MPI_Send( &data[nucs_starting_pos], last_nucs, MPI_CHAR,
           an_id, send_data_tag, MPI_COMM_WORLD);


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
       }
       else
       {
         //nucs_starting_pos = reads_this_stage;
         //reads_starting_pos = read_count;
         printf("reads this stage = %d     current step = %d    for id = %d\n", reads_this_stage, read_count, an_id);
         printf("nucleotide starting position = %d     reads starting position = %d\n", nucs_starting_pos, reads_starting_pos);
         //send info message
         ierr = MPI_Send( &read_count, 1 , MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);

         ierr = MPI_Send( &data_pos[reads_starting_pos], read_count, MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);
         printf("FINISHED SENDING POSITIONAL READ DATA\n");

         //SEND THE STRING OF READS
         ierr = MPI_Send( &reads_this_stage, 1 , MPI_INT,
           an_id, send_data_tag, MPI_COMM_WORLD);

         ierr = MPI_Send( &data[nucs_starting_pos], reads_this_stage, MPI_CHAR,
           an_id, send_data_tag, MPI_COMM_WORLD);


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


         read_count += avg_reads_per_process;
         an_id++;
         nucs_starting_pos = reads_this_stage;
         reads_starting_pos += avg_reads_per_process;
       }

    }

  }
  */      

  }

//##########################################################################
//STUFF SENT TO OTHER PROCESSORS  
  else
  {
    //printf("RECEIVING DATA\n");
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    //RECEIVE DATA FOR READS and READ START AND STOPPING LOCATIONS
    ierr = MPI_Recv( &avg_reads_to_receive, 1, MPI_INT, 
               root_process, send_data_tag, MPI_COMM_WORLD, &status);
    ierr = MPI_Recv( &data_pos_local, avg_reads_to_receive, MPI_INT, 
               root_process, send_data_tag, MPI_COMM_WORLD, &status);
    //printf("RECEIVING CONCATENATED READ DATA\n");
    ierr = MPI_Recv( &reads_to_receive, 1, MPI_INT,
               root_process, send_data_tag, MPI_COMM_WORLD, &status);

    ierr = MPI_Recv( &data_local, reads_to_receive, MPI_CHAR,
    	    root_process, send_data_tag, MPI_COMM_WORLD, &status);
    
    //printf("RECEIVED CONCATENATED READ DATA\n");
    // RECEIVE DATA FOR MSA
    ierr = MPI_Recv( &avg_seq_len_local, 1, MPI_INT,
               root_process, send_data_tag, MPI_COMM_WORLD, &status);

    ierr = MPI_Recv( &total_seq_count_local, 1, MPI_INT,
               root_process, send_data_tag, MPI_COMM_WORLD, &status);
    
    //printf("RECEIVED SEQUENCE INFORMATION\n");
//###################################################################################
//RECEIVE MSA INDIVIDUAL SEQUENCES AND BEGIN LOOPING READS OVER THEM
    //printf("BEGINNING TO RECEIVE MSA SEQUENCES AND PERFORMING READ MAPPING\n");    
    for(n = 0; n < total_seq_count_local; n++)
    {
	    
      ierr = MPI_Recv( &msa_data_local, avg_seq_len_local, MPI_CHAR,
   	      root_process, send_data_tag, MPI_COMM_WORLD, &status);
       
      char *str = msa_data_local;
      int seq_size = strlen(str);

      /* LOOP THROUGH READS AND BEGIN CUTTING INTO KMERS */
      printf("BEGIN LOOPING THROUGH READS WITH PROC %d\n", my_id);
      int pos = 0;
      int counted_nucs = 0;
      for( z = 0; z < avg_reads_to_receive - 2; z++)
      {
	loop = 0;
	//printf("%d\n", data_pos_local[z]);
	pos = data_pos_local[z];
	//########################################################
        for( k = 0; k < pos; k++)
        {
          read_buffer[k] = data_local[counted_nucs + k];
          
	  //printf("%c", read_buffer[k]);
        }
	counted_nucs += pos;
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
	    //printf("%c", ps);
        
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
	  //printf("\n");
          read_hash_array[0][j] = read_hash_number;
          read_hash_number = 0;
        }


        //printf("***SEQUENCE HASHING SLIDING WINDOW***\n");
        int msa_hash_number = 0;
	//printf("msa_hash_number = %d\n", msa_hash_number);
        for(i = 0 ; i < seq_size; i++)
        {
          //printf("%c\n", str[i]);
          int msa_hash = 0;
	  //hash = 0;
	  //printf("reguar msa hash = %d\n", msa_hash);
	  switch(str[i])
          {
            case 'A':
            case 'a':
              msa_hash = 0;
              break;
            case 'C':
            case 'c':
              msa_hash = 1;
              break;
            case 'G':
            case 'g':
              msa_hash = 2;
              break;
            case 'T':
            case 't':
              msa_hash = 3;
              break;
	    default :
	      msa_hash = 4;
          }
	  //printf("reguar msa hash = %d    nuc = %c", msa_hash, str[i]);
	  int sub_hash = 0;
	  //printf("SUB HASH = %d\n", sub_hash);
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
	  //printf("sub_hash = %d\n", sub_hash);
          msa_hash_number = msa_hash_number * 4 + msa_hash - (sub_hash << 28);
          for(j = 0; j < read_km_number; j++)
          {
            //printf("read kmer = %d   msa kmer = %d\n", read_hash_array[0][j], msa_hash_number);
            if(read_hash_array[0][j] == msa_hash_number)
            {
//############################################################
              printf("FOUND HASH MATCH seq position:%d in loci: %d read kmer:%d of read %d     at ID = %d\n", i, n, j, z, my_id);
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
	//printf("msa hash = %d\n", hash);
	msa_hash = 0;
	sub_hash = 0;
	  
//####################################################################

//####################################################################
      }
      msa_hash_number = 0;
//FINISH LOOPING OVER MSA AND DOING COMPARISON
//#####################################################################

//#####################################################################
//FINISH LOOPING OVER CONTIGUOUS READ STRING AND MOVE TO NEXT READ
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

//END OF INDIVIDUAL PROCESSOR/CORE JOBS
//############################################################################

  

ierr = MPI_Finalize();

  return 0;
}




