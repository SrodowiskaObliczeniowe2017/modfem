/************************************************************************
File pcs_mpi_intf.c
Contains definitions of routines: 
  pcr_init_parallel - to initialize parallel communication library
  pcr_my_proc_id - to get id of current process
  pcr_nr_proc - to get number of all processes  

  pcr_send_buffer_open - to open a buffer for sending data
  pcr_buffer_pack_int - to pack an array of integers to a buffer
  pcr_buffer_pack_double - to pack an array of doubles to a buffer
  pcr_buffer_pack_char - to pack an array of characters to a buffer
  pcr_buffer_send - to send a buffer to destination processor
  !pcr_buffer_start_send - to start sending a buffer to destination processor
  !pcr_buffer_finish_send - to finish sending a buffer to destination processor
  pcr_buffer_bcast - to broadcast the content of a buffer 

  pcr_buffer_receive - to receive a buffer 
  !pcr_buffer_start_receive - to start receiving a buffer
  !pcr_buffer_finish_receive - to finish receiving a buffer 
  pcr_buffer_unpack_int - to unpack an array of integers from a buffer
  pcr_buffer_unpack_double - to unpack an array of doubles from a buffer
  pcr_buffer_unpack_char - to unpack an array of characters from a buffer
  pcr_recv_buffer_close - to close a buffer for receiveing data

  pcr_send_int - to send without buffering an array of integers
  pcr_send_double - to send without buffering an array of doubles 
  pcr_receive_int - to receive without buffering an array of integers
  pcr_receive_double - to receive without buffering an array of doubles 
  pcr_bcast_double - to broadcast an array of doubles
  pcr_bcast_int - to broadcast an array of integers
  pcr_bcast_char - to broadcast an array of characyers
  pcr_allreduce_sum_int - to reduce (sum) and broadcast an array of integers 
  pcr_allreduce_sum_double - to reduce (sum) and broadcast an array of doubles
  pcr_allreduce_max_int - to reduce (max) and broadcast an array of integers 
  pcr_allreduce_max_double - to reduce (max) and broadcast an array of doubles
  pcr_exit_parallel - to exit parallel communication library
  pcr_barrier - to synchronise paralell execution
  pcr_is_parallel_initialized - to determine wheter parallelism is initialized or not

REMARK: the module is prepared to serve multiple buffers; 

PROCESSORS (PROCESSES) ARE NUMBERED FROM 1 to NR_PROC !!!!!!!!!!!!

------------------------------  			
History:        
	02.2002 - Krzysztof Banas, initial version	
	11.2005 - Krzysztof Koszyka, version 1.1
    07.2009 - Krzysztof Banas, version 2.0
    10.2012 - Kazimierz Michalik, version 2.1
*************************************************************************/
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include<assert.h>
#include<stdint.h>

/* parallel communication interface specification */
#include "../include/pch_intf.h"
#include "../include/uth_log.h" //#include "../include/dbg.h"

/*** TYPEDEFS ***/
typedef struct {
  int buffer_id;
  int message_id;
  int buffer_size;
  void *buffer;
  int sender_id; // in MPI style -> 0 <= sender_id < nr_proc
  int receiver_id; // in MPI style -> 0 <= sender_id < nr_proc
  int position;
  int data_size;
  MPI_Request request;
  MPI_Status status;
} pct_buffer_struct;

/*** CONSTANTS ***/
const int PCC_ANY_PROC = -1; /* wildcard for arbitrary processor id */
const int PCC_USE_CURRENT_BUFFER = -1;
const int PCC_MASTER_PROC_ID = 1;

#define MAX_NR_BUFFERS  1

FILE *output_stream;

const int PCC_DEFAULT_BUFFER_SIZE=1000000; /* 1 MB, 125 000 doubles */

// for implementation with buffered sends an internal MPI buffer is created
void *MPI_buffer;
// it is assumed that at most ten LARGE outgoing messages are buffered
// i.e. each subdomain exchanges data with maximum 10 neighbouring subdomains
int MPI_buffer_size;

pct_buffer_struct buffer_array[MAX_NR_BUFFERS];
int pcv_nr_buffers = 0;

// my_proc_id = my_rank + 1 !!!!!!
int pcv_my_rank;

/*** GLOBAL VARIABLES ***/
int pcv_nr_proc;
int pcv_my_proc_id;
 

//------------------------------------------------------//
// pcr_is_parallel_initialized - to determine wheter parallelism is initialized or not
//------------------------------------------------------//
extern int pcr_is_parallel_initialized( //returns 1 - is initialized, 0 - is NOT initialized
void
)
{
  int is_mpi_initialized=0;
  return ( MPI_SUCCESS == MPI_Initialized(&is_mpi_initialized) );
}

//------------------------------------------------------//
// pcr_my_proc_rank - to get rank of current process (id = rank +1)
//------------------------------------------------------//
int pcr_my_proc_rank( //returns id of current process
 							  void
 						   )
{
   return pcv_my_rank;
}
 
//------------------------------------------------------//
// pcr_my_proc_id - to get id of current process (id = rank +1)
//------------------------------------------------------//
int pcr_my_proc_id( //returns id of current process
 							  void
 						   )
{
   return pcv_my_proc_id;
}
 
//------------------------------------------------------//
// pcr_nr_proc - to get number of all processes  
//------------------------------------------------------//
int pcr_nr_proc( //returns number of all process
 					  void
 					   )
{
  return pcv_nr_proc;
}
 


/*---------------------------------------------------------
  pcr_init_parallel - to initialize parallel communication library
---------------------------------------------------------*/
int pcr_init_parallel( /* returns: >=0 - success code, <0 - error code */
  int* argc,       /* in: */
  char** argv,     /* in: */
  char* Work_dir,  /* in: */
  char *interactive_output_name, /* in/out: */
  FILE **interactive_output_p, /* out: reset pointer to interactive output */
  int* Nr_pr,      /* out: number of processors/processes */
  int* My_id       /* out: local process ID */
  )
{
  
  pct_buffer_struct* buffer_p;
  int i, rc;
  FILE *interactive_output;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  MPI_Init(argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&pcv_my_rank);
  /*!!!!!!!!!!!!!!!!*/
  pcv_my_proc_id = pcv_my_rank+1;
  /*!!!!!!!!!!!!!!!!*/
  MPI_Comm_size(MPI_COMM_WORLD,&pcv_nr_proc);
  
/*kbw*/
  printf("Initialized MPI_Comm: number of process(es)/(ors) %d, my id %d\n",
  pcv_nr_proc,pcv_my_proc_id);
/*kew*/

/*kbw
  printf("output file: %s\n", interactive_output_name);
/*kew*/

  // set interactive output for parallel execution
  if(*interactive_output_p != stdout){
/*kbw
  printf("output file: %s\n", interactive_output_name);
/*kew*/
    sprintf(interactive_output_name,"%s_%d",
        interactive_output_name,pcv_my_proc_id);
/*kbw
  printf("output file: %s\n", interactive_output_name);
/*kew*/
    interactive_output = fopen(interactive_output_name,"w");
  }
  else interactive_output = *interactive_output_p;

  output_stream = interactive_output;

  if(interactive_output == NULL){
    fprintf(stderr, "Cannot establish interactive output. Exiting\n");
    exit(0);
  }

#ifdef DEBUG_PCM
  fprintf(interactive_output,
	  "Starting parallel communication module in debug mode.\n");
#endif

  
  *Nr_pr = pcv_nr_proc;
  *My_id = pcv_my_proc_id;

  pcv_my_proc_id = pcv_my_proc_id;
  pcv_nr_proc = pcv_nr_proc;
  
  for(i=0;i<MAX_NR_BUFFERS;i++){
    buffer_array[i].buffer_id = i; 
    buffer_array[i].message_id = -1;
    buffer_array[i].buffer = NULL;
    buffer_array[i].buffer_size = 0;
  }

  // initialize internal MPI buffer
  // it is assumed that at most ten LARGE outgoing messages are buffered
  // i.e. each subdomain exchanges data with maximum 10 neighbouring subdomains
  MPI_buffer_size = 10*PCC_DEFAULT_BUFFER_SIZE;
  MPI_buffer = (void*)malloc(MPI_buffer_size);
  rc = MPI_Buffer_attach(MPI_buffer, MPI_buffer_size);
  //log_info("MPI buffer size set to %d Bytes",MPI_buffer_size);
  if (rc != MPI_SUCCESS) {
    fprintf(interactive_output,
	    "Buffer attach failed. Return code= %d Terminating\n", rc);
    MPI_Finalize();
    exit(0);
  }
  
  *interactive_output_p = interactive_output;
  
  return(0);
}


/*---------------------------------------------------------
pcr_print_master - to return the id of read/write master processor
---------------------------------------------------------*/
int pcr_print_master(
  /* returns the id of read/write master processor or error code <0 */
){

  return(PCC_MASTER_PROC_ID);
}


/*---------------------------------------------------------
  pcr_send_buffer_open - to open a buffer for sending data
---------------------------------------------------------*/
int pcr_send_buffer_open( 
  /* returns: >=0 - message buffer ID = success code, <0 - error code */
  int Message_id,  /* in: message ID */
  int Buffer_size  /* in: message buffer size (0 for DEFAULT_BUFFER_SIZE) */
  )
{

  pct_buffer_struct* buffer_p;
  int i;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* check the number of opened buffers */
  pcv_nr_buffers++;

#ifdef DEBUG_PCM
  if(pcv_nr_buffers > MAX_NR_BUFFERS){
    fprintf(output_stream, "Too much MPI buffers requested. \n");
    fprintf(output_stream, "Increase MAX_NR_BUFFERS in pcs_mpi_intf.c.\n");
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif  

  /* set current_buffer to the first free buffer */
  for(i=0;i<MAX_NR_BUFFERS;i++){
    if(buffer_array[i].message_id == -1){
      buffer_p = &buffer_array[i];
      buffer_p->buffer_id = i;
      break;
    }
  }
  
  buffer_p->message_id = Message_id;
  
  /* initialize communication buffer data structure */
  
  if(Buffer_size == 0) Buffer_size = PCC_DEFAULT_BUFFER_SIZE;

  if(buffer_p->buffer_size < Buffer_size){
    if(buffer_p->buffer == NULL){
      buffer_p->buffer_size = Buffer_size;
      buffer_p->buffer = (void*) malloc( (unsigned)buffer_p->buffer_size );
    } 
    else {
      buffer_p->buffer_size = Buffer_size;
      buffer_p->buffer = realloc( buffer_p->buffer,  
				  (unsigned)buffer_p->buffer_size);
      fprintf(output_stream, "Warning: reallocating buffer in MPI!\n");
    }
  } 

  buffer_p->sender_id = pcv_my_proc_id-1;
  buffer_p->receiver_id = -1;
  buffer_p->message_id = Message_id;
  buffer_p->position = 0;
#ifdef DEBUG_PCM
  buffer_p->data_size = 0;
#endif
  
  return(buffer_p->buffer_id);
}


/*---------------------------------------------------------
  pcr_buffer_pack_int - to pack an array of integers to a buffer
---------------------------------------------------------*/
int pcr_buffer_pack_int( /* returns: >=0 - success code, <0 - error code */
  const int Message_id,  /* in: message ID */
  const int Buffer_id,  /* in: buffer ID */
  const int Nr_num,     /* in: number of values to pack to the buffer */
  const int* Numbers    /* in: array of numbers to pack */
)
{
  int rozmiar=0;
  pct_buffer_struct* buffer_p;
  int i;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* pointer to buffer data structure */
  buffer_p = &buffer_array[Buffer_id];

#ifdef DEBUG_PCM  
  if(Message_id != buffer_p->message_id){
    fprintf(output_stream, "Error - packing to different message !\n");
    fprintf(output_stream, "buffer message_id %d, input message__id %d\n", 
	    buffer_p->message_id,Message_id);
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	

#ifdef DEBUG_PCM  
  MPI_Pack_size(Nr_num, MPI_INT, MPI_COMM_WORLD, &rozmiar);
  buffer_p->data_size+=rozmiar;
  
  if(buffer_p->data_size > buffer_p->buffer_size) {
    fprintf(output_stream, "Too large message %d. Increase buffer_size!\n",
	    Message_id);
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	
  
  MPI_Pack((void*)Numbers, Nr_num, MPI_INT, buffer_p->buffer, 
	   buffer_p->buffer_size, 
           &buffer_p->position, MPI_COMM_WORLD);

#ifdef DEBUG_PCM  
  /* check integrity of data */
  if(buffer_p->position != buffer_p->data_size){
    fprintf(output_stream, "Error in packing MPI buffer in pcr_buffer_pack_double!\n");
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	

  return(0);
}

/*---------------------------------------------------------
  pcr_buffer_pack_double - to pack an array of doubles to a buffer
---------------------------------------------------------*/
int pcr_buffer_pack_double( /* returns: >=0 - success code, <0 - error code */
  const int Message_id,  /* in: message ID */
  const int Buffer_id,  /* in: message ID */
  const int Nr_num,     /* in: number of values to pack to the buffer */
  const double* Numbers /* in: array of numbers to pack */
)
{
  int rozmiar=0;
  pct_buffer_struct* buffer_p;
  int i;

/*++++++++++++++++ executable statements ++++++++++++++++*/


/* pointer to buffer data structure */
//  assert(Buffer_id < pcv_nr_buffers);
  buffer_p = &buffer_array[pcv_nr_buffers-1];
  
#ifdef DEBUG_PCM  
  if(Message_id != buffer_p->message_id){
    fprintf(output_stream, "Error - packing to different message !\n");
    fprintf(output_stream, "buffer message_id %d, input message__id %d\n", 
	    buffer_p->message_id,Message_id);
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	

#ifdef DEBUG_PCM  
  MPI_Pack_size(Nr_num, MPI_DOUBLE, MPI_COMM_WORLD, &rozmiar);
  buffer_p->data_size+=rozmiar;
  
  if(buffer_p->data_size > buffer_p->buffer_size) {
    fprintf(output_stream, "Too large message %d. Increase buffer_size!\n",
	    Message_id);
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	
  
  MPI_Pack((void*)Numbers, Nr_num, MPI_DOUBLE, buffer_p->buffer, 
	   buffer_p->buffer_size, 
           &buffer_p->position, MPI_COMM_WORLD);

#ifdef DEBUG_PCM  
  /* check integrity of data */
  if(buffer_p->position != buffer_p->data_size){
    fprintf(output_stream, "Error in packing MPI buffer in pcr_buffer_pack_double!\n");
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	

  return(0);
}

/*---------------------------------------------------------
  pcr_buffer_pack_char - to pack an array of chars to a buffer
---------------------------------------------------------*/
int pcr_buffer_pack_char( /* returns: >=0 - success code, <0 - error code */
  const int Message_id,  /* in: message ID */
  const int Buffer_id,  /* in: buffer ID */
  const int Nr_num,     /* in: number of values to pack to the buffer */
  const char* Numbers /* in: array of numbers to pack */
)
{
  int rozmiar=0;
  pct_buffer_struct* buffer_p;
  int i;

/*++++++++++++++++ executable statements ++++++++++++++++*/


/* pointer to buffer data structure */
  buffer_p = &buffer_array[pcv_nr_buffers-1];
  
#ifdef DEBUG_PCM  
  if(Message_id != buffer_p->message_id){
    fprintf(output_stream, "Error - packing to different message !\n");
    fprintf(output_stream, "buffer message_id %d, input message__id %d\n", 
	    buffer_p->message_id,Message_id);
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	

#ifdef DEBUG_PCM  
  MPI_Pack_size(Nr_num, MPI_INT, MPI_COMM_WORLD, &rozmiar);
  buffer_p->data_size+=rozmiar;
  
  if(buffer_p->data_size > buffer_p->buffer_size) {
    fprintf(output_stream, "Too large message %d. Increase buffer_size!\n",
	    Message_id);
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	
  
  MPI_Pack((void*)Numbers, Nr_num, MPI_CHAR, buffer_p->buffer, 
	   buffer_p->buffer_size, 
           &buffer_p->position, MPI_COMM_WORLD);

#ifdef DEBUG_PCM  
  /* check integrity of data */
  if(buffer_p->position != buffer_p->data_size){
    fprintf(output_stream, "Error in packing MPI buffer in pcr_buffer_pack_char!\n");
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	

  return(0);
}


/*---------------------------------------------------------
  pcr_buffer_send - to send a buffer to destination processor
---------------------------------------------------------*/
int pcr_buffer_send( /* returns: >=0 - success code, <0 - error code */
  int Message_id,     /* in: message ID */
  int Buffer_id,  /* in: buffer ID */
  int Dest_proc   /* in : destination processor ID */
  )
{   

  pct_buffer_struct* buffer_p;
  MPI_Request *request; 
/*++++++++++++++++ executable statements ++++++++++++++++*/


/* pointer to buffer data structure */
  buffer_p = &buffer_array[Buffer_id];

#ifdef DEBUG_PCM  
  if(Message_id != buffer_p->message_id){
    fprintf(output_stream, "Error - sending different message !\n");
    fprintf(output_stream, "buffer message_id %d, input message__id %d\n", 
	    buffer_p->message_id,Message_id);
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	

  buffer_p->receiver_id = Dest_proc-1;

/*kbw
  fprintf(output_stream, "Sending message %d to proc %d\n", 
	  buffer_p->message_id, buffer_p->receiver_id);
/*kew*/

  /* buffered blocking send */
  MPI_Bsend(buffer_p->buffer, buffer_p->position, MPI_PACKED,
  	   buffer_p->receiver_id, buffer_p->message_id, MPI_COMM_WORLD);
  //MPI_Isend(buffer_p->buffer, buffer_p->position, MPI_PACKED,
  //	    buffer_p->receiver_id, buffer_p->message_id, MPI_COMM_WORLD,request);

  /* prepare buffer for new sends or receives */
  pcv_nr_buffers--;

  buffer_p->message_id = -1;
  buffer_p->sender_id = -1;
  buffer_p->receiver_id = -1;
  buffer_p->position = 0;
#ifdef DEBUG_PCM  
  buffer_p->data_size = 0;
#endif

  return(0);
}


/*--------------------------------------------------------
  pcr_buffer_bcast - to broadcast a buffer to all processes
---------------------------------------------------------*/
int pcr_buffer_bcast( /* returns: >=0 - success code, <0 - error code*/
  int Message_id, /* in: message ID */
  int Buffer_id,  /* in: buffer ID */
  int Sender_proc /* in: sender processor ID */
  )    
{

  pct_buffer_struct* buffer_p;
  int i;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  
  if(pcv_my_proc_id==Sender_proc){
    /* pointer to buffer data structure */
    buffer_p = &buffer_array[Buffer_id];
  }

#ifdef DEBUG_PCM  
  if(pcv_my_proc_id==Sender_proc){
  if(Message_id != buffer_p->message_id){
    fprintf(output_stream, "Error - broadcasting different message !\n");
    fprintf(output_stream, "buffer message_id %d, input message__id %d\n", 
	    buffer_p->message_id,Message_id);
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
  }
#endif	

#ifdef DEBUG_PCM  
  if(pcv_my_proc_id==Sender_proc){
    /* check integrity of data */
    if(buffer_p->position != buffer_p->data_size){
      fprintf(output_stream, "Error in packing MPI buffer in pcr_buffer_bcast!\n");
      fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
      exit(0);
    }
  }
#endif	


  if(pcv_my_proc_id!=Sender_proc){

    /* check the number of opened buffers */
    pcv_nr_buffers++;
    
#ifdef DEBUG_PCM
    if(pcv_nr_buffers > MAX_NR_BUFFERS){
      fprintf(output_stream, "Too much MPI buffers requested. \n");
      fprintf(output_stream, "Increase MAX_NR_BUFFERS in pcs_mpi_intf.c.\n");
      fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
      exit(0);
    }
#endif
    
    /* set current_buffer to the first free buffer */
    for(i=0;i<MAX_NR_BUFFERS;i++){
      if(buffer_array[i].message_id == 0){
	buffer_p = &buffer_array[i];
	buffer_p->buffer_id = i;
	break;
      }
    }
    
    buffer_p->message_id = Message_id;
    
  }

  MPI_Bcast(&buffer_p->position, 1, MPI_INT, Sender_proc-1, MPI_COMM_WORLD);

  if(pcv_my_proc_id!=Sender_proc){

    /* initialize communication buffer data structure */

    if(buffer_p->buffer_size < buffer_p->position) 
      buffer_p->buffer_size = buffer_p->position;  
    if(buffer_p->buffer_size < PCC_DEFAULT_BUFFER_SIZE)
      buffer_p->buffer_size = PCC_DEFAULT_BUFFER_SIZE;

    if(buffer_p->buffer == NULL){
      buffer_p->buffer = (void*) malloc( (unsigned)buffer_p->buffer_size );
    } 
    else {
      buffer_p->buffer = realloc( buffer_p->buffer,  
				  (unsigned)buffer_p->buffer_size);
      fprintf(output_stream, "Warning: reallocating buffer in MPI!\n");
    }
    
    
    /* initialize communication buffer data structure */
    buffer_p->message_id = Message_id;
    buffer_p->sender_id = Sender_proc-1;
    buffer_p->receiver_id = pcv_my_proc_id-1;
    buffer_p->position = 0;
#ifdef DEBUG_PCM
    buffer_p->data_size = 0;
#endif
    
  }	

  MPI_Bcast(buffer_p->buffer, buffer_p->position, MPI_PACKED,
	      buffer_p->sender_id, MPI_COMM_WORLD);
   
  
  if(pcv_my_proc_id==Sender_proc){

    /* prepare buffer for new sends or receives */
    pcv_nr_buffers--;
    
    buffer_p->sender_id = -1;
    buffer_p->receiver_id = -1;
    buffer_p->message_id = -1;
    buffer_p->position = 0;
#ifdef DEBUG_PCM
    buffer_p->data_size = 0;
#endif
    
  } else {
    
    /* prepare buffer for unpacking */
    buffer_p->position = 0;
    
  }

  return (0);
}


/*---------------------------------------------------------
  pcr_buffer_receive - to receive a buffer from a particular processor
---------------------------------------------------------*/
int pcr_buffer_receive( 
		   /* returns: >= 0 - Buffer_id, <0 - error code */
  int Message_id,   /* in: id of the message containing the data */
  int Sender_proc,  /* in : sender processor ID */
  int Buffer_size /* in: message buffer size (0 for DEFAULT_BUFFER_SIZE) */
  )
{

  pct_buffer_struct* buffer_p;
  int i;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* check the number of opened buffers */
  pcv_nr_buffers++;

#ifdef DEBUG_PCM
  if(pcv_nr_buffers > MAX_NR_BUFFERS){
    fprintf(output_stream, "Too much MPI buffers requested. \n");
    fprintf(output_stream, "Increase MAX_NR_BUFFERS in pcs_mpi_intf.c.\n");
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif
  
  /* set current_buffer to the first free buffer */
  for(i=0;i<MAX_NR_BUFFERS;i++){
    if(buffer_array[i].message_id == -1){
      buffer_p = &buffer_array[i];
      buffer_p->buffer_id = i;
      break;
    }
  }
  
  buffer_p->message_id = Message_id;
  
/*kbw
  fprintf(output_stream, "found buffer %d to receive message %d\n",
	  buffer_p->buffer_id, buffer_p->message_id);
/*kew*/


  /* initialize communication buffer data structure */
  buffer_p->message_id = Message_id;
  if(Sender_proc == PCC_ANY_PROC) buffer_p->sender_id = MPI_ANY_SOURCE;
  else buffer_p->sender_id = Sender_proc-1;
  buffer_p->receiver_id = pcv_my_proc_id-1;
  buffer_p->position = 0;
#ifdef DEBUG_PCM
  buffer_p->data_size = 0;
#endif

/*kbw
  fprintf(output_stream, "receiving buffer parameters:\n");
  fprintf(output_stream, "ID %d, message ID %d, sender_id %d, receiver_id %d\n",
	  buffer_p->buffer_id, buffer_p->message_id,
	  buffer_p->sender_id, buffer_p->receiver_id);
/*kew*/


  if(Buffer_size == 0) Buffer_size = PCC_DEFAULT_BUFFER_SIZE;

  int result = MPI_Probe(buffer_p->sender_id, buffer_p->message_id,  
	    MPI_COMM_WORLD, &buffer_p->status);

  if(result!=MPI_SUCCESS) fprintf(output_stream, "Error in MPI_Probe!!!!!!\n");

  MPI_Get_count(&buffer_p->status, MPI_PACKED, &i);

  if(Buffer_size < i) Buffer_size = i;

  if(buffer_p->buffer_size < Buffer_size){
    buffer_p->buffer_size = Buffer_size;
    if(buffer_p->buffer == NULL){
      buffer_p->buffer = (void*) malloc( (unsigned)buffer_p->buffer_size );
    } 
    else {
      buffer_p->buffer = realloc( buffer_p->buffer,  
				  (unsigned)buffer_p->buffer_size);
      fprintf(output_stream, "Warning: reallocating buffer in MPI!\n");
    }
  } 
  
/*kbw
  fprintf(output_stream, "after MPI_Probe:\n");
  fprintf(output_stream, "ID %d, message ID %d, size %d\n",
	  buffer_p->buffer_id, buffer_p->message_id,
	  buffer_p->buffer_size);
/*kew*/

#ifdef DEBUG_PCM
  buffer_p->data_size = buffer_p->buffer_size;
#endif

  MPI_Recv(buffer_p->buffer, buffer_p->buffer_size, MPI_PACKED,
	   buffer_p->sender_id, buffer_p->message_id, MPI_COMM_WORLD, 
	   &buffer_p->status);

/*kbw
  fprintf(output_stream, "after MPI_Recv:\n");
/*kew*/

  return(buffer_p->buffer_id);
}


/*---------------------------------------------------------
  pcr_buffer_unpack_int - to unpack an array of integers from a buffer
---------------------------------------------------------*/
int pcr_buffer_unpack_int( /* returns: >=0 - success code, <0 - error code */
  int Message_id,  /* in: message ID */
  int Buffer_id,  /* in: buffer ID */
  int Nr_num,     /* in: number of values to pack to the buffer */
  int* Numbers    /* in: array of numbers to pack */
)
{

  pct_buffer_struct* buffer_p;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* pointer to buffer data structure */
  buffer_p = &buffer_array[Buffer_id];

#ifdef DEBUG_PCM  
  if(Message_id != buffer_p->message_id){
    fprintf(output_stream, "Error - unpacking from different message !\n");
    fprintf(output_stream, "buffer message_id %d, input message__id %d\n", 
	    buffer_p->message_id,Message_id);
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	

  MPI_Unpack(buffer_p->buffer, buffer_p->buffer_size, &buffer_p->position,
	     Numbers, Nr_num, MPI_INT, MPI_COMM_WORLD);

  return(0);
}

/*---------------------------------------------------------
  pcr_buffer_unpack_double - to unpack an array of integers from a buffer
---------------------------------------------------------*/
int pcr_buffer_unpack_double( /* returns: >=0 - success code, <0 - error code */
  int Message_id,  /* in: message ID */
  int Buffer_id,  /* in: buffer ID */
  int Nr_num,     /* in: number of values to unpack from the buffer */
  double* Numbers /* in: array of numbers to unpack */
)
{

  pct_buffer_struct* buffer_p;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* pointer to buffer data structure */
  buffer_p = &buffer_array[Buffer_id];

#ifdef DEBUG_PCM  
  if(Message_id != buffer_p->message_id){
    fprintf(output_stream, "Error - unpacking from different message !\n");
    fprintf(output_stream, "buffer message_id %d, input message__id %d\n", 
	    buffer_p->message_id,Message_id);
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	


  MPI_Unpack(buffer_p->buffer, buffer_p->buffer_size, &buffer_p->position,
	     Numbers, Nr_num, MPI_DOUBLE, MPI_COMM_WORLD);
	
  return(0);
}

/*---------------------------------------------------------
  pcr_buffer_unpack_char - to unpack an array of integers from a buffer
---------------------------------------------------------*/
int pcr_buffer_unpack_char( /* returns: >=0 - success code, <0 - error code */
  int Message_id,  /* in: message ID */
  int Buffer_id,  /* in: buffer ID */
  int Nr_num,     /* in: number of values to unpack from the buffer */
  char* Numbers /* in: array of numbers to unpack */
)
{

  pct_buffer_struct* buffer_p;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* pointer to buffer data structure */
  buffer_p = &buffer_array[Buffer_id];

#ifdef DEBUG_PCM  
  if(Message_id != buffer_p->message_id){
    fprintf(output_stream, "Error - unpacking from different message !\n");
    fprintf(output_stream, "buffer message_id %d, input message__id %d\n", 
	    buffer_p->message_id,Message_id);
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	


  MPI_Unpack(buffer_p->buffer, buffer_p->buffer_size, &buffer_p->position,
	     Numbers, Nr_num, MPI_CHAR, MPI_COMM_WORLD);
	
  return(0);
}

/*---------------------------------------------------------
  pcr_recv_buffer_close - to close a buffer for receiveing data
---------------------------------------------------------*/
int pcr_recv_buffer_close( /* returns: >=0 - success code, <0 - error code */
  int Message_id,  /* in: message ID */
  int Buffer_id  /* in: buffer ID */
  )
{

  pct_buffer_struct* buffer_p;

/*++++++++++++++++ executable statements ++++++++++++++++*/
  
#ifdef DEBUG_PCM
/* check the number of opened buffers */
  if(pcv_nr_buffers == 0){
    fprintf(output_stream, "Closing non existing communication buffer. \n");
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif

/* pointer to buffer data structure */
  buffer_p = &buffer_array[Buffer_id];

#ifdef DEBUG_PCM  
  if(Message_id != buffer_p->message_id){
    fprintf(output_stream, "Error - packing to different message !\n");
    fprintf(output_stream, "buffer message_id %d, input message__id %d\n", 
	    buffer_p->message_id,Message_id);
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif	


  /* clear communication buffer data structure */
  buffer_p->message_id = -1;
  buffer_p->buffer_id = -1;
  buffer_p->sender_id = -1;
  buffer_p->receiver_id = -1;
  buffer_p->position = 0;
  buffer_p->data_size = 0;

  pcv_nr_buffers--;

  return(0);

}

/*---------------------------------------------------------
  pcr_send_int - to send an array of integers
---------------------------------------------------------*/
extern int pcr_send_int( 
                    /* returns: >=0 - success code, <0 - error code */
  const int Dest_proc_id, /* in : destination processor ID */
  const int Message_id,   /* in: id of the message containing the data */
  const int Nr_num,       /* in: number of values to pack to the buffer */
  const int* Numbers      /* in: array of numbers to pack */
)
{
  //log_info("pcr_send_int: sending %d ints from %d to %d",Nr_num,pcr_my_proc_id(),Dest_proc_id);
  assert(Numbers != NULL);
  return MPI_Send((void*)Numbers, Nr_num, MPI_INT,
		  Dest_proc_id-1, Message_id, MPI_COMM_WORLD);
}

/*---------------------------------------------------------
  pcr_send_double - to send an array of doubles 
---------------------------------------------------------*/
extern int pcr_send_double( 
                    /* returns: >=0 - success code, <0 - error code */
  const int Dest_proc_id, /* in : destination processor ID */
  const int Message_id,   /* in: id of the message containing the data */
  const int Nr_num,       /* in: number of values to send */
  const double* Numbers   /* in: array of numbers to send */
  )
{
  assert(Numbers != NULL);
  return MPI_Send((void*)Numbers, Nr_num, MPI_DOUBLE,
		  Dest_proc_id-1, Message_id, MPI_COMM_WORLD);
}

/*---------------------------------------------------------
  pcr_receive_int - to receive an array of integers
---------------------------------------------------------*/
extern int pcr_receive_int( 
                    /* returns: >=0 - success code, <0 - error code */
  int Sender_proc_id, /* in : destination processor ID */
  int Message_id,   /* in: id of the message containing the data */
  int Nr_num,       /* in: number of values to pack to the buffer */
  int* Numbers      /* in: array of numbers to pack */
  )
{
  assert(Numbers != NULL);
  //log_info("pcr_receive_int: reciving %d ints from %d to %d",Nr_num,Sender_proc_id,pcr_my_proc_id());
  MPI_Status status;
  return MPI_Recv((void*)Numbers, Nr_num, MPI_INT, Sender_proc_id-1, Message_id,
		  MPI_COMM_WORLD, &status);
}

/*---------------------------------------------------------
  pcr_receive_double - to receive an array of doubles 
---------------------------------------------------------*/
extern int pcr_receive_double( 
                    /* returns: >=0 - success code, <0 - error code */
  int Sender_proc_id, /* in : destination processor ID */
  int Message_id,   /* in: id of the message containing the data */
  int Nr_num,       /* in: number of values to pack to the buffer */
  double* Numbers   /* in: array of numbers to pack */
  )
{
  assert(Numbers != NULL);
  MPI_Status status;
  return MPI_Recv((void*)Numbers, Nr_num, MPI_DOUBLE, Sender_proc_id-1, Message_id,
		  MPI_COMM_WORLD, &status);
}

/*---------------------------------------------------------
  pcr_send_long - to send an array of integers
---------------------------------------------------------*/
extern int pcr_send_long(
                    /* returns: >=0 - success code, <0 - error code */
  int Dest_proc_id, /* in : destination processor ID */
  int Message_id,   /* in: id of the message containing the data */
  int Nr_num,       /* in: number of values to send */
  long int* Numbers   /* in: array of numbers to send */
  )
{
  assert(Numbers != NULL);
  return MPI_Send((void*)Numbers, Nr_num, MPI_LONG,
		  Dest_proc_id-1, Message_id, MPI_COMM_WORLD);
}

/*---------------------------------------------------------
  pcr_receive_long - to receive an array of integers
---------------------------------------------------------*/
extern int pcr_receive_long(
                    /* returns: >=0 - success code, <0 - error code */
  int Sender_proc_id, /* in : destination processor ID */
  int Message_id,   /* in: id of the message containing the data */
  int Nr_num,       /* in: number of values to pack to the buffer */
  long int* Numbers      /* in: array of numbers to pack */
  )
{
  assert(Numbers != NULL);
  MPI_Status status;
  return MPI_Recv((void*)Numbers, Nr_num, MPI_LONG, Sender_proc_id-1, Message_id,
		  MPI_COMM_WORLD, &status);
}

/*---------------------------------------------------------
  pcr_send_bytes - to send an array of integers
---------------------------------------------------------*/
extern int pcr_send_bytes(
                    /* returns: >=0 - success code, <0 - error code */
  const int Dest_proc_id, /* in : destination processor ID */
  const int Message_id,   /* in: id of the message containing the data */
  const int Nr_bytes,       /* in: number of values to send */
  const uint8_t* Bytes   /* in: array of numbers to send */
  )
{
  assert(Bytes != NULL);
  return MPI_Send((void*)Bytes, Nr_bytes, MPI_BYTE,
		  Dest_proc_id-1, Message_id, MPI_COMM_WORLD);
}

/*---------------------------------------------------------
  pcr_receive_bytes - to receive an array of integers
---------------------------------------------------------*/
extern int pcr_receive_bytes(
                    /* returns: >=0 - success code, <0 - error code */
  int Sender_proc_id, /* in : destination processor ID */
  int Message_id,   /* in: id of the message containing the data */
  int Nr_bytes,       /* in: number of values to pack to the buffer */
  uint8_t* Bytes      /* in: array of numbers to pack */
  )
{
  assert(Bytes != NULL);
  MPI_Status status;
  return MPI_Recv((void*)Bytes, Nr_bytes, MPI_BYTE, Sender_proc_id-1, Message_id,
		  MPI_COMM_WORLD, &status);
}


/*--------------------------------------------------------
  pcr_bcast_double - to broadcast an array of doubles
---------------------------------------------------------*/
int pcr_bcast_double( /* returns: >=0 - success code, <0 - error code*/
  int Sender_proc_id, /* in : source processor ID */
  int Nr_num,       /* in: number of values to pack to the buffer */
  double* Numbers   /* in: array of numbers to pack */
  )    
{

  MPI_Bcast((void*)Numbers,Nr_num,MPI_DOUBLE,Sender_proc_id-1,MPI_COMM_WORLD);
  
  return (0);
}


/*--------------------------------------------------------
  pcr_bcast_int - to broadcast an array of integers
---------------------------------------------------------*/
int pcr_bcast_int( /* returns: >=0 - success code, <0 - error code*/
  int Sender_proc_id, /* in : source processor ID */
  int Nr_num,       /* in: number of values to pack to the buffer */
  int* Numbers      /* in: array of numbers to pack */
  )    
{

  MPI_Bcast((void*)Numbers,Nr_num,MPI_INT,Sender_proc_id-1,MPI_COMM_WORLD);

  return (0);
}


/*--------------------------------------------------------
  pcr_bcast_char - to broadcast an array of charegers
---------------------------------------------------------*/
int pcr_bcast_char( /* returns: >=0 - success code, <0 - error code*/
  int Sender_proc_id, /* in : source processor ID */
  int Nr_num,       /* in: number of values to pack to the buffer */
  char* Numbers      /* in: array of numbers to pack */
  )    
{

  MPI_Bcast((void*)Numbers,Nr_num,MPI_CHAR,Sender_proc_id-1,MPI_COMM_WORLD);

  return (0);
}


/*---------------------------------------------------------
pcr_allreduce_sum_int - to reduce by summing and broadcast an array of integers 
---------------------------------------------------------*/
int pcr_allreduce_sum_int( /* returns: >=0 - success code, <0 - error code*/
  const int Nr_num,          /* in: number of values to reduce */
  const int* Numbers,        /* in: array of numbers to sum */
   int* Numbers_reduced /* in: array of summed numbers */
  )    
{

  assert(Numbers != Numbers_reduced);
  //log_info("pcr_allreduce_sum_int");
  MPI_Allreduce((void*)Numbers,(void*)Numbers_reduced,Nr_num,
		MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        
  return (0);
}


/*--------------------------------------------------------
pcr_allreduce_sum_double - to reduce by summing and broadcast an array of doubles
---------------------------------------------------------*/
int pcr_allreduce_sum_double(/* returns: >=0 - success code, <0 - error code*/
  const int Nr_num,          /* in: number of values to reduce */
  const double* Numbers,        /* in: array of numbers to sum */
  double *Numbers_reduced /* in: array of summed numbers */
  )
{
  assert(Numbers != Numbers_reduced);
  //log_info("pcr_allreduce_sum_double");
  MPI_Allreduce((void*)Numbers,(void*)Numbers_reduced,Nr_num,
		MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        
 return (0);
}

/*---------------------------------------------------------
pcr_allreduce_max_int - to reduce by maxming and broadcast an array of integers 
---------------------------------------------------------*/
int pcr_allreduce_max_int( /* returns: >=0 - success code, <0 - error code*/
  const int Nr_num,          /* in: number of values to reduce */
  const int* Numbers,        /* in: array of numbers to max */
  int* Numbers_reduced /* in: array of maxmed numbers */
  )    
{
  assert(Numbers != Numbers_reduced);
  //log_info("pcr_allreduce_max_int");
  MPI_Allreduce((void*)Numbers,(void*)Numbers_reduced,Nr_num,
		MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  
  return (0);
}


/*--------------------------------------------------------
pcr_allreduce_max_double - to reduce by maxming and broadcast an array of doubles
---------------------------------------------------------*/
int pcr_allreduce_max_double( /* returns: >=0 - success code, <0 - error code*/
  int Nr_num,          /* in: number of values to reduce */
  double* Numbers,        /* in: array of numbers to max */
  double* Numbers_reduced /* in: array of maxmed numbers */
  )    
{
  assert(Numbers != Numbers_reduced);
  //log_info("pcr_allreduce_max_double");
  MPI_Allreduce((void*)Numbers,(void*)Numbers_reduced,Nr_num,
		MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        
 return (0);
}

/*---------------------------------------------------------
  pcr_exit_parallel - to exit parallel communication library
---------------------------------------------------------*/
int pcr_exit_parallel( /* returns: >=0 - success code, <0 - error code */
  )
{
    
  pct_buffer_struct* buffer_p;
  int i;

#ifdef DEBUG_PCM
  if(pcv_nr_buffers > 0){
    fprintf(output_stream, "Leaving MPI with opened communication buffers.\n");
    fprintf(output_stream, "My_proc_id %d. Exiting.\n", pcv_my_proc_id);
    exit(0);
  }
#endif

  for(i=0;i<MAX_NR_BUFFERS;i++){
    if(buffer_array[i].buffer != NULL) free(buffer_array[i].buffer);
  }

  MPI_Buffer_detach(&MPI_buffer, &MPI_buffer_size);
  free(MPI_buffer);

/*kbw
  printf("Exiting MPI: %d\n");
/*kew*/

  return ( MPI_Finalize() );
}

//------------------------------------------------------//
// pcr_barrier - to synchronise parallel execution
//------------------------------------------------------//
int pcr_barrier( // returns >=0 - success, <0- error code
						 void
						 )
{
  return MPI_Barrier(MPI_COMM_WORLD);
}

//------------------------------------------------------//
// pcr_allgather_int - gathers together values from a group of processors and distributes to all
//------------------------------------------------------//
int pcr_allgather_int( //returns success code
                              const int send_values[],
                              const int n_send_values,
                              int gathered_values[],
                              int n_gathered_values

)
{
    //log_info("pcr_allgather_int");
    return MPI_Allgather((void*)send_values,n_send_values,MPI_INT,
                         (void*)gathered_values,n_gathered_values,MPI_INT,
			 MPI_COMM_WORLD);
}


