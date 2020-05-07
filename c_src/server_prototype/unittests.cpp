#include <mpi.h>
#include <string>
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <stdlib.h>
#include <vector>
#include <stdarg.h>
#include "serv_util.h"

using ::testing::Each;
using ::testing::Eq;

//must be called by all ranks in world in the same file and line with the same expression
//tests if the expression is true
#define SYNCH_TEST(expr) synchTest(expr, __FILE__, __LINE__, #expr)

//no restriction on when and where to call it: ranks do not need to synchronize
//tests if the expression is true
#define LONE_TEST(expr) loneTest(expr, __FILE__, __LINE__, #expr)

//called by several (not necessarily all) ranks to ensure they have the same value of something
//tests if the expression evaluates to the same long int on all participating ranks.
#define RANK_COMPARE_TEST(expr) rankCompareTest(expr, __FILE__, __LINE__, #expr)

#define ABORT_ON_NULL(ptr, message) if(ptr==NULL){nullAbort(__FILE__, __LINE__, #ptr, message);}else{}


//describes a position in the source code, with a filename, linenumber and expression (as a name)
class Position {

public:

	char file_name[100];
	char var_name[100];
	int line_no;

	bool operator <(const Position &other) const{
		int file_comp = strcmp(file_name, other.file_name);
		if (file_comp != 0){
			//files are different, order by file
			return file_comp < 0;
		}
		else{
			//Positions in same file, order by line number
			if (line_no != other.line_no){
				//different lines, order by that
				return line_no < other.line_no;
			}
			else{
				//Same file, same line, go by variable name
				int var_comp = strcmp(var_name, other.var_name);
				return var_comp < 0;
			}
		}
	}

	friend std::ostream & operator<< (std::ostream &out, const Position &pos){

		out << "\"" << pos.var_name << "\"" << " in line "
			<< pos.line_no << " of file " << pos.file_name;

		return out;

	}

};

//invalid before PositionedValue_mpi_create()
MPI_Datatype Position_mpi;

void Position_mpi_create(){

	int blockcounts[3] = {100,100,1};
	MPI_Datatype types[3] = {MPI_CHAR, MPI_CHAR, MPI_INT};
	MPI_Aint displs[3];
	
	struct Position p;
	
	MPI_Get_address(&p.file_name, &displs[0]);
	MPI_Get_address(&p.var_name, &displs[1]);
	MPI_Get_address(&p.line_no, &displs[2]);

	size_t start_addr = displs[0];
	for (int i=0; i<3; i++){
		displs[i] -= start_addr;
	}

	MPI_Type_create_struct(3,blockcounts, displs, types, &Position_mpi);
	MPI_Type_commit(&Position_mpi);

}

// prints a vector of tuples, where the first item represents a rank
// and the second a value taken at that rank.
// Do this in table format, so that it is clear which value is taken at which rank
std::ostream& operator<<(std::ostream& out, std::vector<std::pair<int,long>> const& values)
{
	out << "Worldrank: Value" << std::endl;
	for (unsigned int i=0; i<values.size(); i++) {
		out << values[i].first << ": " <<values[i].second << std::endl;
	}
	return out;
}

// Test fixture to run a test with several MPI-ranks, enables the use of the x_Test macros defined above
class CollectiveTest : public testing::Test {

	int comm_rank, comm_size;

	std::vector<int> synch_test_results;
	std::vector<std::string> synch_test_names;
	
	std::string lone_test_report;

	std::vector<Position> rank_compare_positions;
	std::vector<long> rank_compare_values;

public:
	CollectiveTest(){
		MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
		MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	}

	int getCommRank(){
		return comm_rank;
	}

	int getCommSize(){
		return comm_size;
	}

	//meant to be used through the SYNCH_TEST macro
	void synchTest(int isPassing, const char *file, int line_no, const char *expr){
		synch_test_results.push_back(isPassing);
		char name[500];
		sprintf(name, "\"%s\" in line %d of file %s", expr, line_no, file);
		synch_test_names.push_back(name);
	}

	//meant to be used through the LONE_TEST macro
	bool loneTest(int isPassing, const char *file, int line_no, const char *expr){
		if (!isPassing){
			char name[500];
			sprintf(name, "\"%s\" failed in line %d of file %s. (rank %d)\n", expr, line_no, file, comm_rank);
			lone_test_report.append(name);
		}
		return isPassing;
	}

	//meant to be used through the RANK_COMPARE_TEST macro
	void rankCompareTest(long value_on_rank, const char *file, int line_no, const char *expr){

		ASSERT_LT(strlen(file), 100) << "File names longer than a 100 characters not supported";
		Position *pos;
		pos = (Position *) calloc(sizeof(Position), 1);
		strcpy(pos->file_name, file);
		strcpy(pos->var_name, expr);
		pos->line_no = line_no;

		rank_compare_positions.push_back(*pos);
		rank_compare_values.push_back(value_on_rank);

	}

protected:

	void TearDown() override{

		checkLoneTests();
		checkSynchTests();
		checkCompareTests();

		MPI_Barrier(MPI_COMM_WORLD);

	}

	
		
private:

	void checkLoneTests(){

		int *report_lengths = NULL;
		if ( comm_rank == 0 ) {
			report_lengths = (int *) malloc(sizeof(int) * comm_size);
			ABORT_ON_NULL(report_lengths, "Malloc Error");
		}
		
		unsigned int my_report_len = lone_test_report.size();
		MPI_Gather(&my_report_len , 1, MPI_UNSIGNED,
				report_lengths, 1, MPI_UNSIGNED,
				0, MPI_COMM_WORLD);

		int *dspl = NULL;
		char *all_reports = NULL;
		if ( comm_rank == 0 ){

			dspl = (int *)malloc(comm_size*sizeof(int));
			ABORT_ON_NULL(dspl, "Mallo Error");
			int total = 0;

			for (int r=0; r<comm_size; r++){
				dspl[r] = total;
				total += report_lengths[r];
				total++; //make room for the null-termination
			}	
			all_reports = (char *)calloc(total,sizeof(char));
		}

		MPI_Gatherv(lone_test_report.c_str(), my_report_len, MPI_CHAR,
				all_reports, report_lengths, dspl, MPI_CHAR,
				0, MPI_COMM_WORLD);
		
		if (comm_rank == 0){
			//iterate over ranks
			for (int r=0; r<comm_size; r++){
				// empty report means no failures. Otherwise, display report
				EXPECT_EQ(0,report_lengths[r]) << "-----Failure in rank " << 
					r << ":\n" << &all_reports[dspl[r]] << std::endl;
			}
		}


		free(all_reports);
		free(report_lengths);
	}


	void checkSynchTests(){
	
		bool synched = synchTestUsedCorrectly();
		if (comm_rank==0){
			ASSERT_TRUE(synched) << "All ranks must call synchTest the same amount of times";
		}

		int *recv = NULL;
		if (comm_rank==0){
			recv = (int *) malloc(synch_test_results.size() * sizeof(int) * comm_size);
			ABORT_ON_NULL(recv, "Mallo Error");
		}

		unsigned long num_tests = synch_test_results.size();
			
		MPI_Gather(&synch_test_results[0], num_tests, MPI_INT,
				recv, num_tests, MPI_INT,
				0, MPI_COMM_WORLD);

		if (comm_rank==0){

			// iterate over all synchtests
			for (unsigned int t=0; t< num_tests; t++){
				std::vector<int> t_fails;
				//iterate over ranks
				for (int r=0; r<comm_size; r++){
					if (!recv[r*num_tests + t]){
						//rank r failed test t
						t_fails.push_back(r);
					}

				}
				char err_mess[500];
				makeErrMsg(err_mess, t_fails, synch_test_names[t].c_str());
				
				EXPECT_EQ(0, t_fails.size()) << err_mess;
			}
		}

		free(recv);
	}

	void checkCompareTests(){


		int *size_on_rank = NULL;
		if (comm_rank==0){
			size_on_rank = (int *) malloc(sizeof(int) * comm_size);
			ABORT_ON_NULL(size_on_rank, "Malloc Error");
		}
		int my_size = rank_compare_values.size();
		MPI_Gather(&my_size, 1, MPI_INT,
				size_on_rank, 1, MPI_INT,
				0, MPI_COMM_WORLD);

		int *displs = NULL;
		long *all_values = NULL;
		Position *all_positions = NULL;

		if (comm_rank==0){
			displs = (int *) malloc(sizeof(int) * comm_size);
			ABORT_ON_NULL(displs, "Malloc Error");
			int total = 0;
			for (int r=0; r<comm_size; r++){
				displs[r] = total;
				total += size_on_rank[r];
			}
			all_values= (long *) malloc(sizeof(long) * total);
			ABORT_ON_NULL(all_values, "Malloc Error");
			all_positions = (Position *) malloc(sizeof(Position) * total);
			ABORT_ON_NULL(all_positions, "Malloc Error");
		}


		// collect values
		long * values = rank_compare_values.size() > 0 ? &rank_compare_values[0] : (long*)MPI_BOTTOM;
		MPI_Gatherv(values, rank_compare_values.size(), MPI_LONG,
				all_values, size_on_rank, displs, MPI_LONG,
				0, MPI_COMM_WORLD);

	
		// collect Positions
		Position* pos = rank_compare_positions.size() > 0 ? &rank_compare_positions[0] : (Position*)MPI_BOTTOM;
		MPI_Gatherv(pos , rank_compare_positions.size(), Position_mpi,
				all_positions, size_on_rank, displs, Position_mpi,
				0, MPI_COMM_WORLD);

	

		if (comm_rank==0){
			// (filename, line_num, expr) -> [rank, value]
			std::map<Position, std::vector<std::pair<int,long>>> map;

			// iterate over ranks r and values v
			for (int r=0; r<comm_size; r++){
				for (int v=0; v<size_on_rank[r]; v++){

					auto pos = all_positions[displs[r]+v];
					auto value = all_values[displs[r]+v];
					std::pair<int, long> rank_value = {r, value};
					auto p = map.find(pos);

					if ( p == map.end()){
						// create new entry for this position
						std::vector<std::pair<int,long>> vec {rank_value};
						map.insert({pos, vec});
					}
					else{
						//add to vector in this position
						p->second.push_back(rank_value);
					}
				}
			}


			// for all (pos, list) entries of the map if the whole list has the same value
			for (auto m_it = map.begin(); m_it != map.end(); m_it++){

				auto vec = m_it->second;

				std::vector<long> plain_values;
				for (unsigned int i=0; i < vec.size(); i++){
					plain_values.push_back(vec[i].second);
				}
				EXPECT_THAT(plain_values, Each(Eq(vec[0].second)))
					<< std::endl
					<< m_it->first << std::endl
					<< "expect all ranks to have the same value, but values are:"
					<< vec;
			}
		}

		free(all_values);
		free(all_positions);
		free(size_on_rank);
		free(displs);

	}
	
	//checks if all ranks called synchTest the same number of times,
	//return is only meaningful in the root process, but must be called by all processes
	bool synchTestUsedCorrectly(){

		int msg_size = synch_test_results.size();

		int *size_on_rank = NULL;
		if (comm_rank==0){
			size_on_rank = (int *) malloc(comm_size * sizeof(int));
			ABORT_ON_NULL(size_on_rank, "Mallo Error");
		}

		MPI_Gather(&msg_size,1,MPI_INT,
				size_on_rank,1,MPI_INT,
				0,MPI_COMM_WORLD);
		if (comm_rank==0){

			int root_size = msg_size;
			for (int i=1; i<comm_size; i++){
				if (root_size != size_on_rank[i]){
					// rank i has different size
					free(size_on_rank);
					return false;
				}
			}
		}	
		free(size_on_rank);
		return true;
	}

	// takes a (possibly empty) vector of ints that indicate on which ranks a test failed
	// and generates an appropriate error message
	void makeErrMsg(char * err_mess, std::vector<int>& t_fails, const char* name){

		if (t_fails.size() == (unsigned long) comm_size){
			sprintf(err_mess, "Test ( %s ) failed on all ranks", name);
		}
		else if (t_fails.size() == 1){
			sprintf(err_mess,
					"Test ( %s ) failed on rank %d",
					name, t_fails[0]);

		}
		else if (t_fails.size() > 0) {
			sprintf(err_mess,
					"Test ( %s ) failed on %ld ranks, the first one being %d",
					name, t_fails.size(), t_fails[0]);
		}
		else {
			sprintf(err_mess,
					"Test ( %s ) succeeded on all ranks",
					name);
		}
	}
};

//////////// TESTS

#include "server.h"

TEST_F(CollectiveTest, t2create_1_server){

	rank_info info;
	int err;
	err = create_servers(MPI_COMM_WORLD, 1, &info);

	SYNCH_TEST(err==0);

	bool is_server = (info.sim_server_rank == 0);

	if (is_server){
		LONE_TEST(info.my_type_size == 1);
		LONE_TEST(info.my_type_rank == 0);
	}
	else{
		LONE_TEST(info.my_type_size == getCommSize()-1);
		RANK_COMPARE_TEST(info.n_sim_ranks);
	}
}

TEST_F(CollectiveTest, create_2_server){

	if (getCommSize() < 4){
		// todo: this can be done more nicely with GTEST_SKIP in
		// gtest version > 1.10
		if (getCommRank() == 0){
			std::cout << "SKIPPED create_2_servers which requires at least 4 ranks\n";
		}
		return ;
	}
	
	rank_info info;
	int err;
	err = create_servers(MPI_COMM_WORLD, 2, &info);

	SYNCH_TEST(err==0);

	bool is_server = (info.sim_server_rank == 0);

	//missing_test: list-test the rank numbers

	if (is_server){
		LONE_TEST(info.my_type_size == 2);
		LONE_TEST(info.my_type_rank < info.my_type_size);
		LONE_TEST(info.my_type_rank >= 0);
	}
	else{
		LONE_TEST(info.my_type_size == getCommSize()-2);
		LONE_TEST(info.my_type_rank < info.my_type_size);
		LONE_TEST(info.my_type_rank >= 0);
	}
}

TEST_F(CollectiveTest, t2create_too_many_servers){

	int num_serv = getCommSize()/2 + 1;
	int err;
	rank_info info;
	err = create_servers(MPI_COMM_WORLD, num_serv, &info);

	SYNCH_TEST(err!=0);
}

TEST_F(CollectiveTest, t2cant_start_simrank_as_server){

	int err;
	rank_info info;
	err = create_servers(MPI_COMM_WORLD, 1, &info);

	LONE_TEST(err == 0);

	bool isServer = (0==info.sim_server_rank);

	if (!isServer){
		struct tasks t;
		int err = run_server(&t, 5, &info, printf);
		LONE_TEST(err!=0);
	}
}

struct PrintfString{
	
	static char* out;

	static void init(int size){
		maxlen = size;
		if (out != NULL){
			throw std::runtime_error("PrintfString still in use, cannot init");
		}
		out = (char*)malloc(sizeof(char)*maxlen);
		out[0] = '\0';
		ABORT_ON_NULL(out, "Malloc error");
		offset = 0;
		
	}

	static int print(const char * format, ...){
		va_list args;
		va_start(args, format);
		int written = vsprintf(out+offset, format, args);
		if (written < 0){
			throw std::runtime_error("vsprintf failed");
		}
		offset += written;
		if (offset >= maxlen){
			throw std::runtime_error("string too small");
		}
		va_end(args);
		return written;
	}

	static void clear(){
		offset = 0;
		free(out);
		out = NULL;
		maxlen = 0;
	}

private:
	static int offset;
	static int maxlen;
};

char* PrintfString::out = NULL;
int   PrintfString::offset = 0;
int   PrintfString::maxlen = 0;

void mock_sim_rank(int *x, int n, MPI_Comm sim_server_comm){

	//sends n integers in x to the server
	//once
	
	MPI_Request reqs[2];

	// tell server about size
	MPI_Igather(&n, 1, MPI_INT,
			NULL, 0, MPI_INT,
			0, sim_server_comm, &reqs[0]);

	//send server the actual data
	MPI_Igatherv(x, n, MPI_INT,
			NULL, NULL, NULL, MPI_INT,
			0, sim_server_comm, &reqs[1]);

	// don't care about the status, but the request must be completed
	// (MPI does not allow mixing synchronous and asynchronous collectives)
	MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);
}

class ServerTest : public CollectiveTest{


	void SetUp() override{
		int err;
		err = create_servers(MPI_COMM_WORLD, 1, &info);
	

		SYNCH_TEST(err==0);

	}

public:
	rank_info info;

};

class TwoServerTest : public CollectiveTest{

	void SetUp() override{
		int err;
		err = create_servers(MPI_COMM_WORLD, 2, &info);

		SYNCH_TEST(err==0);

	}

public:
	rank_info info;


};

class MassServerTest : public CollectiveTest{

	void SetUp() override{

		int err;
		//make half the ranks into servers
		int n_servs = getCommSize()/2;
		err = create_servers(MPI_COMM_WORLD, n_servs, &info);

		SYNCH_TEST(err==0);

	}

public:
	rank_info info;


};



TEST_F(ServerTest, t2simple_server_calc){

	//only average, but every step
	struct tasks tsk = {1,0,0};
	
	bool isServer = (0==info.sim_server_rank);

	if (isServer){
		PrintfString::init(10000);
		// one timestep, only analyze
		run_server(&tsk, 1, &info, PrintfString::print);
		char expected[50];
		sprintf(expected, "0: mean=%lf\n", (info.n_sim_ranks+1)/2.0);
		if ( ! LONE_TEST(strcmp(PrintfString::out, expected)==0 ) ){
			std::cout << "server output wrong. Was:\n"
				<<PrintfString::out << std::endl;
			std::cout << "server output should be:\n"
				<<expected << std::endl;
		}
		PrintfString::clear();
	}
	else{
		int m = info.sim_server_rank;
		int x[3] = {m, m, m};
		mock_sim_rank(x, 3, info.sim_server_comm);
	}
}

TEST_F(ServerTest, t2uneven_sim_extreme){

	//only average, but every step
	struct tasks tsk = {1,0,0};
	
	bool isServer = (0==info.sim_server_rank);

	if (isServer){
		PrintfString::init(10000);
		// one timestep, only analyze
		run_server(&tsk, 1, &info, PrintfString::print);
		char expected[50];
		sprintf(expected, "0: mean=%lf\n", (5+1)/2.0);
		if ( ! LONE_TEST(strcmp(PrintfString::out, expected)==0 ) ){
			std::cout << "server output wrong. Was:\n"
				<<PrintfString::out << std::endl;
			std::cout << "server output should be:\n"
				<<expected << std::endl;
		}
		PrintfString::clear();
	}
	else{
		if (info.sim_server_rank == 1){
			int x[5];
			for (int i=0; i<5; i++){
				x[i] = i+1;
			}
			mock_sim_rank(x, 5, info.sim_server_comm);
		}
		else{
			mock_sim_rank(NULL, 0, info.sim_server_comm);
		}
	}
}

TEST_F(ServerTest, partly_uneven_sim){

	//only average, but every step
	struct tasks tsk = {1,0,0};
	
	bool isServer = (0==info.sim_server_rank);

	if (info.n_sim_ranks < 2) {
		if (getCommRank() == 0){
			std::cout << "not enough simranks, skipping test" << std::endl;
		}
		return;

	}

	if (isServer){
		PrintfString::init(10000);
		// one timestep, only analyze
		run_server(&tsk, 1, &info, PrintfString::print);
		char expected[50];
		sprintf(expected, "0: mean=%lf\n", (10+1)/2.0);
		if ( ! LONE_TEST(strcmp(PrintfString::out, expected)==0 ) ){
			std::cout << "server output wrong. Was:\n"
				<<PrintfString::out << std::endl;
			std::cout << "server output should be:\n"
				<<expected << std::endl;
		}
		PrintfString::clear();
	}
	else{
		if (info.sim_server_rank == 1){
			int x[7];
			for (int i=0; i<7; i++){
				x[i] = i+1;
			}
			mock_sim_rank(x, 7,  info.sim_server_comm);
		}
		else if (info.sim_server_rank == 2){
			int x[3];
			for (int i=0; i<3; i++){
				x[i] = i+7+1;
			}
			mock_sim_rank(x, 3,  info.sim_server_comm);
		}
		else{
			mock_sim_rank(NULL, 0,  info.sim_server_comm);
		}
	}
}

TEST_F(ServerTest, t2multi_time_step){

	//only calc on every other timestep
	struct tasks tsk = {2,0,0};

	bool isServer = (info.sim_server_rank == 0);

	int size0 = 100;
	int size2 = 50;
	int size4 = 200;

	if (isServer){
		PrintfString::init(10000);
		// one timestep, only analyze
		run_server(&tsk, 5, &info, PrintfString::print);
		char expected[50];
		double exval0 = (size0+1)/2.0;
		double exval2 = (size2+1)/2.0;
 		double exval4 = (size4+1)/2.0;
		sprintf(expected, "0: mean=%lf\n2: mean=%lf\n4: mean=%lf\n",
				exval0, exval2, exval4);
		if ( ! LONE_TEST(strcmp(PrintfString::out, expected)==0 ) ){
			std::cout << "server output wrong. Was:\n"
				<<PrintfString::out << std::endl;
			std::cout << "server output should be:\n"
				<<expected << std::endl;
		}
		PrintfString::clear();
	}
	else{
		std::vector<int> sizes = {size0, size2, size4};
		for (int size : sizes){
			int *x = (int *)malloc(sizeof(int)*size);
			ABORT_ON_NULL(x, "Malloc Error");
			for (int i=0; i<size; i++){
				x[i] = i+1;
			}
			int size_per_rank = size/info.my_type_size + 1;
			int start = std::min(size_per_rank*(info.my_type_rank), size);
			int end = std::min(size_per_rank*(info.my_type_rank+1), size);
			mock_sim_rank(x+start, end-start, info.sim_server_comm);
			free(x);
		}
	}
}

TEST_F(ServerTest, t2variance0){
	struct tasks tsk = {0,1,0};

	bool isServer = (info.sim_server_rank == 0);

	if (isServer){
		PrintfString::init(10000);
		run_server(&tsk, 1, &info, PrintfString::print);
		char expected[50];
		sprintf(expected, "0: var=0.000000\n");

		if ( ! LONE_TEST(strcmp(PrintfString::out, expected)==0)){
			std::cout << "server output wrong. Was:\n"
				<<PrintfString::out << std::endl;
			std::cout << "server output should be:\n"
				<<expected << std::endl;
		}
		PrintfString::clear();
		
	}
	else{
		const int size = 100;
		int x[size];
		for (int i=0; i<size; i++){
			x[i] = 7;
		}
		mock_sim_rank(x, size, info.sim_server_comm);
	}
}

TEST_F(ServerTest, t2variance_nonzero){

	struct tasks tsk = {0,1,0};
	
	bool isServer = (info.sim_server_rank == 0);

	if (isServer){
		PrintfString::init(10000);
		run_server(&tsk, 1, &info, PrintfString::print);
		char expected[50];
		sprintf(expected, "0: var=1.000000\n");

		if ( ! LONE_TEST(strcmp(PrintfString::out, expected)==0)){
			std::cout << "server output wrong. Was:\n"
				<<PrintfString::out << std::endl;
			std::cout << "server output should be:\n"
				<<expected << std::endl;
		}
		PrintfString::clear();
	}
	else{
		const int size = 10;
		int x[size];
		for (int i=0; i < size; i+=2){
			x[i] = 1;
			x[i+1] = -1;
		}
		mock_sim_rank(x, size, info.sim_server_comm);
	}
}

TEST_F(TwoServerTest, mean_with_all_ranks_equal){

	struct tasks tsk = {1,0,0};

	const int size = 30;

	bool isServer = (info.sim_server_rank == 0);

	if (isServer){
		PrintfString::init(10000);
		run_server(&tsk, 1, &info, PrintfString::print);	

		char expected[50];

		if (info.my_type_rank == 0){

			sprintf(expected, "0: mean=%lf\n", (size+1)/2.0);
		}
		else{
			//other servers shouldnt print
			sprintf(expected, "\0");

		}
		if ( ! LONE_TEST(strcmp(PrintfString::out, expected)==0)){

			std::cout << "server output wrong. Was:\n"
				<<PrintfString::out << std::endl;
			std::cout << "server output should be:\n"
				<<expected << std::endl;
		}
		PrintfString::clear();

	}
	else{
		int x[size];
		for (int i=0; i<size; i++){
			x[i] = size-i;
		}
		mock_sim_rank(x, size, info.sim_server_comm);
	}
}

TEST_F(TwoServerTest, mean_ranks_different_all_obs){

	const int size = 1017;

	struct tasks tsk = {1,1,1};

	bool isServer = (info.sim_server_rank == 0);

	if (isServer){

		PrintfString::init(10000);
		run_server(&tsk, 1, &info, PrintfString::print);	

		char expected[50];

		if (info.my_type_rank == 0){

			double mean = (size+1)/2.0;
			double var = (size+1)*(2*size+1)/6.0  -  mean*mean;
			double absmax = (double) size;
			sprintf(expected, "0: mean=%lfvar=%lfabsmax=%lf\n",
					mean, var, absmax);
		}
		else{
			//other servers shouldnt print
			sprintf(expected, "\0");

		}
		if ( ! LONE_TEST(strcmp(PrintfString::out, expected)==0)){

			std::cout << "server output wrong. Was:\n"
				<<PrintfString::out << std::endl;
			std::cout << "server output should be:\n"
				<<expected << std::endl;
		}
		PrintfString::clear();

	}
	else{
		int x[size];
		for (int i=0; i<size; i++){
			x[i] = size-i;
		}

		int size_per_rank = size/info.my_type_size + 1;
		int start = std::min(size_per_rank*(info.my_type_rank), size);
		int end = std::min(size_per_rank*(info.my_type_rank+1), size);

		mock_sim_rank(x+start, end-start, info.sim_server_comm);
	}



}

TEST_F(MassServerTest, mean_ranks_different){

	const int size = 1017;

	struct tasks tsk = {1,0,0};

	bool isServer = (info.sim_server_rank == 0);

	if (isServer){

		PrintfString::init(10000);
		run_server(&tsk, 1, &info, PrintfString::print);	

		char expected[50];

		if (info.my_type_rank == 0){

			sprintf(expected, "0: mean=%lf\n", (size+1)/2.0);
		}
		else{
			//other servers shouldnt print
			sprintf(expected, "\0");

		}
		if ( ! LONE_TEST(strcmp(PrintfString::out, expected)==0)){

			std::cout << "server output wrong. Was:\n"
				<<PrintfString::out << std::endl;
			std::cout << "server output should be:\n"
				<<expected << std::endl;
		}
		PrintfString::clear();

	}
	else{
		int x[size];
		for (int i=0; i<size; i++){
			x[i] = size-i;
		}

		int size_per_rank = size/info.my_type_size + 1;
		int start = std::min(size_per_rank*(info.my_type_rank), size);
		int end = std::min(size_per_rank*(info.my_type_rank+1), size);

		mock_sim_rank(x+start, end-start, info.sim_server_comm);
	}

}


//////////// TESTS END

int main(int argc, char**argv){

	int failure;

	::testing::InitGoogleTest(&argc, argv);

	MPI_Init(NULL, NULL);

	Position_mpi_create();

	// only rank 0 should print test results
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	::testing::TestEventListeners& listeners = \
		::testing::UnitTest::GetInstance()->listeners();
	if (rank != 0){
		delete listeners.Release(listeners.default_result_printer());
	}
	
	failure = RUN_ALL_TESTS();

	if (failure != 0){
		DPRINT("nonzero exit code: %d", failure);
	}
	else if (failure == 0 && rank == 0){
		int size;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		std::cout << "tests were successful with " << size << " ranks" << std::endl;
	}

	MPI_Finalize();

	return failure;

}
