#include <soma_config.h>
/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.6
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "SOMA"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "SOMA"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION get_soma_version()
#endif

enum enum_pseudo_random_number_generator { pseudo_random_number_generator__NULL = -1, pseudo_random_number_generator_arg_PCG32 = 0, pseudo_random_number_generator_arg_MT, pseudo_random_number_generator_arg_TT800 };
enum enum_move_type { move_type__NULL = -1, move_type_arg_TRIAL = 0, move_type_arg_SMART };
enum enum_iteration_alg { iteration_alg__NULL = -1, iteration_alg_arg_POLYMER = 0, iteration_alg_arg_SET };
enum enum_set_generation_algorithm { set_generation_algorithm__NULL = -1, set_generation_algorithm_arg_SIMPLE = 0, set_generation_algorithm_arg_FIXEDMINUS_NMINUS_SETS };

/** @brief Where the command line options are stored */
struct som_args
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *detailed_help_help; /**< @brief Print help, including all details and hidden options, and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * coord_file_arg;	/**< @brief File containing the system description and all Coordinates. (HDF5-Format).  */
  char * coord_file_orig;	/**< @brief File containing the system description and all Coordinates. (HDF5-Format) original value given at command line.  */
  const char *coord_file_help; /**< @brief File containing the system description and all Coordinates. (HDF5-Format) help description.  */
  int timesteps_arg;	/**< @brief Number of MC sweeps carried out by SOMA..  */
  char * timesteps_orig;	/**< @brief Number of MC sweeps carried out by SOMA. original value given at command line.  */
  const char *timesteps_help; /**< @brief Number of MC sweeps carried out by SOMA. help description.  */
  char * ana_file_arg;	/**< @brief File containing the analysis frequency and is going to be appended by the new measured observables. (HDF5-Format) (default='').  */
  char * ana_file_orig;	/**< @brief File containing the analysis frequency and is going to be appended by the new measured observables. (HDF5-Format) original value given at command line.  */
  const char *ana_file_help; /**< @brief File containing the analysis frequency and is going to be appended by the new measured observables. (HDF5-Format) help description.  */
  int gpus_arg;	/**< @brief Number of GPUs per MPI-node to use. The devices 0-(gpus-1) are going to be occupied by SOMA. Every node must feature this number of nodes and it is assumed that the scheduler assignes the ranks consequently to the nodes. If set to 0, SOMA tries to run on host-code. Ignored if compiled without OPENACC. (default='0').  */
  char * gpus_orig;	/**< @brief Number of GPUs per MPI-node to use. The devices 0-(gpus-1) are going to be occupied by SOMA. Every node must feature this number of nodes and it is assumed that the scheduler assignes the ranks consequently to the nodes. If set to 0, SOMA tries to run on host-code. Ignored if compiled without OPENACC. original value given at command line.  */
  const char *gpus_help; /**< @brief Number of GPUs per MPI-node to use. The devices 0-(gpus-1) are going to be occupied by SOMA. Every node must feature this number of nodes and it is assumed that the scheduler assignes the ranks consequently to the nodes. If set to 0, SOMA tries to run on host-code. Ignored if compiled without OPENACC. help description.  */
  int only_gpu_arg;	/**< @brief Specify a specific Device for all ranks. Useful for MPI-single rank runs. This option overrides the --gpus option..  */
  char * only_gpu_orig;	/**< @brief Specify a specific Device for all ranks. Useful for MPI-single rank runs. This option overrides the --gpus option. original value given at command line.  */
  const char *only_gpu_help; /**< @brief Specify a specific Device for all ranks. Useful for MPI-single rank runs. This option overrides the --gpus option. help description.  */
  double screen_output_interval_arg;	/**< @brief Specify the number of seconds between an output about timings on the screen. (default='10').  */
  char * screen_output_interval_orig;	/**< @brief Specify the number of seconds between an output about timings on the screen. original value given at command line.  */
  const char *screen_output_interval_help; /**< @brief Specify the number of seconds between an output about timings on the screen. help description.  */
  int rng_seed_arg;	/**< @brief Global seed for the pseudo-random number generator. If you pass seed < 0, seed = time(NULL) will be used. Option useful for debuggin purposes. (default='-1').  */
  char * rng_seed_orig;	/**< @brief Global seed for the pseudo-random number generator. If you pass seed < 0, seed = time(NULL) will be used. Option useful for debuggin purposes. original value given at command line.  */
  const char *rng_seed_help; /**< @brief Global seed for the pseudo-random number generator. If you pass seed < 0, seed = time(NULL) will be used. Option useful for debuggin purposes. help description.  */
  enum enum_pseudo_random_number_generator pseudo_random_number_generator_arg;	/**< @brief Option to select the pseudo random number generator. (default='PCG32').  */
  char * pseudo_random_number_generator_orig;	/**< @brief Option to select the pseudo random number generator. original value given at command line.  */
  const char *pseudo_random_number_generator_help; /**< @brief Option to select the pseudo random number generator. help description.  */
  int omp_threads_arg;	/**< @brief Number of omp threads used per MPI rank. If you pass n < 1 it will be set to 1. (default='1').  */
  char * omp_threads_orig;	/**< @brief Number of omp threads used per MPI rank. If you pass n < 1 it will be set to 1. original value given at command line.  */
  const char *omp_threads_help; /**< @brief Number of omp threads used per MPI rank. If you pass n < 1 it will be set to 1. help description.  */
  int nonexact_area51_flag;	/**< @brief Specify to use the exact check of area51. This includes checks, whether a particle moves through a forbidden area51. Performance might be slightly increased if switched to nonexact. Configuration generation is always in exact mode. (default=off).  */
  const char *nonexact_area51_help; /**< @brief Specify to use the exact check of area51. This includes checks, whether a particle moves through a forbidden area51. Performance might be slightly increased if switched to nonexact. Configuration generation is always in exact mode. help description.  */
  enum enum_move_type move_type_arg;	/**< @brief Specify the Monte-Carlo move type. (default='SMART').  */
  char * move_type_orig;	/**< @brief Specify the Monte-Carlo move type. original value given at command line.  */
  const char *move_type_help; /**< @brief Specify the Monte-Carlo move type. help description.  */
  enum enum_iteration_alg iteration_alg_arg;	/**< @brief Specify the iteration algorithm of the beads. This specifies also the level of parallelism that is possible. (default='POLYMER').  */
  char * iteration_alg_orig;	/**< @brief Specify the iteration algorithm of the beads. This specifies also the level of parallelism that is possible. original value given at command line.  */
  const char *iteration_alg_help; /**< @brief Specify the iteration algorithm of the beads. This specifies also the level of parallelism that is possible. help description.  */
  int skip_tests_flag;	/**< @brief Skip tests SOMA is usually preforming before and after the simulation to ensure integrety of the data. (default=off).  */
  const char *skip_tests_help; /**< @brief Skip tests SOMA is usually preforming before and after the simulation to ensure integrety of the data. help description.  */
  int load_balance_arg;	/**< @brief Frequency of the load balancer. For homogenous architectures this can be set to high values, for hetereogenous architectures across the MPI ranks small values help to equilibrate faster. Non-MPI runs are uneffected. Values < 0 deactivate the load-balancer. (default='500').  */
  char * load_balance_orig;	/**< @brief Frequency of the load balancer. For homogenous architectures this can be set to high values, for hetereogenous architectures across the MPI ranks small values help to equilibrate faster. Non-MPI runs are uneffected. Values < 0 deactivate the load-balancer. original value given at command line.  */
  const char *load_balance_help; /**< @brief Frequency of the load balancer. For homogenous architectures this can be set to high values, for hetereogenous architectures across the MPI ranks small values help to equilibrate faster. Non-MPI runs are uneffected. Values < 0 deactivate the load-balancer. help description.  */
  double accepted_load_inbalance_arg;	/**< @brief  [0,100] Percent of step time which is ignored by load balancer. Low values enable better load balancing, but could cause fluctuation of polymers. (default='8').  */
  char * accepted_load_inbalance_orig;	/**< @brief  [0,100] Percent of step time which is ignored by load balancer. Low values enable better load balancing, but could cause fluctuation of polymers. original value given at command line.  */
  const char *accepted_load_inbalance_help; /**< @brief  [0,100] Percent of step time which is ignored by load balancer. Low values enable better load balancing, but could cause fluctuation of polymers. help description.  */
  int autotuner_restart_period_arg;	/**< @brief Period in which the autotuner is restarted. (default='10000').  */
  char * autotuner_restart_period_orig;	/**< @brief Period in which the autotuner is restarted. original value given at command line.  */
  const char *autotuner_restart_period_help; /**< @brief Period in which the autotuner is restarted. help description.  */
  char * user_arg;	/**< @brief Additional arguments. The usage of these arguments defined by the user. The default setting ignores the arguments..  */
  char * user_orig;	/**< @brief Additional arguments. The usage of these arguments defined by the user. The default setting ignores the arguments. original value given at command line.  */
  const char *user_help; /**< @brief Additional arguments. The usage of these arguments defined by the user. The default setting ignores the arguments. help description.  */
  enum enum_set_generation_algorithm set_generation_algorithm_arg;	/**< @brief Option to select the algorithm to generate the indepent sets. (default='SIMPLE').  */
  char * set_generation_algorithm_orig;	/**< @brief Option to select the algorithm to generate the indepent sets. original value given at command line.  */
  const char *set_generation_algorithm_help; /**< @brief Option to select the algorithm to generate the indepent sets. help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int detailed_help_given ;	/**< @brief Whether detailed-help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int coord_file_given ;	/**< @brief Whether coord-file was given.  */
  unsigned int timesteps_given ;	/**< @brief Whether timesteps was given.  */
  unsigned int ana_file_given ;	/**< @brief Whether ana-file was given.  */
  unsigned int gpus_given ;	/**< @brief Whether gpus was given.  */
  unsigned int only_gpu_given ;	/**< @brief Whether only-gpu was given.  */
  unsigned int screen_output_interval_given ;	/**< @brief Whether screen-output-interval was given.  */
  unsigned int rng_seed_given ;	/**< @brief Whether rng-seed was given.  */
  unsigned int pseudo_random_number_generator_given ;	/**< @brief Whether pseudo-random-number-generator was given.  */
  unsigned int omp_threads_given ;	/**< @brief Whether omp-threads was given.  */
  unsigned int nonexact_area51_given ;	/**< @brief Whether nonexact-area51 was given.  */
  unsigned int move_type_given ;	/**< @brief Whether move-type was given.  */
  unsigned int iteration_alg_given ;	/**< @brief Whether iteration-alg was given.  */
  unsigned int skip_tests_given ;	/**< @brief Whether skip-tests was given.  */
  unsigned int load_balance_given ;	/**< @brief Whether load-balance was given.  */
  unsigned int accepted_load_inbalance_given ;	/**< @brief Whether accepted-load-inbalance was given.  */
  unsigned int autotuner_restart_period_given ;	/**< @brief Whether autotuner-restart-period was given.  */
  unsigned int user_given ;	/**< @brief Whether user was given.  */
  unsigned int set_generation_algorithm_given ;	/**< @brief Whether set-generation-algorithm was given.  */

} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure som_args (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure som_args (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *som_args_purpose;
/** @brief the usage string of the program */
extern const char *som_args_usage;
/** @brief the description string of the program */
extern const char *som_args_description;
/** @brief all the lines making the help output */
extern const char *som_args_help[];
/** @brief all the lines making the detailed help output (including hidden options and details) */
extern const char *som_args_detailed_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct som_args *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct som_args *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct som_args *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct som_args *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct som_args *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the detailed help (including hidden options and details)
 */
void cmdline_parser_print_detailed_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed som_args structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct som_args *args_info);
/**
 * Deallocates the string fields of the som_args structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct som_args *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct som_args *args_info,
  const char *prog_name);

extern const char *cmdline_parser_pseudo_random_number_generator_values[];  /**< @brief Possible values for pseudo-random-number-generator. */
extern const char *cmdline_parser_move_type_values[];  /**< @brief Possible values for move-type. */
extern const char *cmdline_parser_iteration_alg_values[];  /**< @brief Possible values for iteration-alg. */
extern const char *cmdline_parser_set_generation_algorithm_values[];  /**< @brief Possible values for set-generation-algorithm. */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
