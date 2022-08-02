ANALYSIS="-*,cert-*,clang-analyzer-*,modernize-*,performance-*,cppcoreguidelines-*,google-*,bugprone-*,misc-*,-google-runtime-references,-cppcoreguidelines-avoid-magic-numbers,-modernize-use-trailing-return-type"

SRC=$(wildcard *.cpp implementation/*.cpp)
OBJ=$(SRC:.cpp=.o)
HDR=$(wildcard *.hh)

STD=c++17
INCLUDES=-I . 
WARNINGS=-Wall -Wextra -Wunused -Wunused-function -Wshadow 
DEBUG=-g -O1

testinterval : $(OBJ)
	$(CXX) $(OBJ) -o $@

%.o : %.cpp $(HDR)
	$(CXX) -MJ $@.json $(INCLUDES) $(WARNINGS) $(DEBUG) -std=$(STD) -c $< -o $@

compile_commands.json :
	echo '[' > compile_commands.json
	cat *.o.json >> compile_commands.json
	echo '{}]' >> compile_commands.json

test: testinterval
	./testinterval

help:
	@echo "make [testinterval] : compile testinterval"
	@echo "make test          : compile and run testinterval"
	@echo "make format        : compile and run testinterval"
	@echo "make analyze       : clang-tidy code analysis"

format :
	clang-format -MJ -i -style=file ${SRC} ${HDR}

clean :
	rm -f testinterval *.o *.o.json *.orig *.plist

analyze :
	clang-tidy --extra-arg="--std=${STD}" -checks=${ANALYSIS} -header-filter=.* ${SRC} -- ${INCLUDES}

fix :
	clang-tidy --fix --extra-arg="--std=${STD}" -checks=${ANALYSIS} -header-filter=.* ${SRC} -- ${INCLUDES}