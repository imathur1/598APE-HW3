FUNC := clang++
copt := -c 
OBJ_DIR := ./bin/
FLAGS := -O3 -march=native -ffast-math -flto -lm -g -Werror -lprofiler -fopenmp
VALGRIND_FLAGS := -O0 -g -fno-omit-frame-pointer -fno-inline -lm -Werror -fopenmp -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free

# Enable profiling via a macro
ifndef PROFILE
	PROFILE := 0
endif
FLAGS += -DPROFILE=$(PROFILE)

CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix $(OBJ_DIR),$(notdir $(CPP_FILES:.cpp=.obj)))

all: main.exe compare_outputs

main.exe:
	$(FUNC) ./main.cpp -o ./main.exe $(FLAGS)

compare_outputs:
	$(FUNC) ./compare_outputs.cpp -o ./compare_outputs $(FLAGS)

test: compare_outputs
	$(FUNC) ./main.cpp -o ./main.exe $(FLAGS) -DPRINT_STEPS=1

clean:
	rm -f ./*.exe ./compare_outputs

view-profile:
	~/go/bin/pprof -http "0.0.0.0:8080" ./main.exe ./my_profile.prof

valgrind: main.cpp
	$(FUNC) ./main.cpp -o ./main.exe $(VALGRIND_FLAGS)

.PHONY: all clean view-profile main.exe compare_outputs test valgrind