CC = g++

EXECS = epr_sim

#Libraries and opt
CFLAGS = -O3 -std=c++11 -static #-fsanitize=address -Wall -O0 -g #DEBUG                  
CFLAGS += -I/opt/local/include  -L/opt/local/lib -lgsl -lgslcblas -lm

all: EPR_sim

EPR_sim: ../src/main_epr.cc
	${CC} -o epr_sim ../src/main_epr.cc ${CFLAGS}

clean:
	rm ${EXECS}
	@echo "Cleaned."
