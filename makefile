BiPoS: Starter.c user_input.c Library.c Synth.c InitFinDistr.c get_params.c Ptr.c Small_routines.c
	gcc -o BiPoS Starter.c user_input.c Library.c Synth.c InitFinDistr.c get_params.c Ptr.c Small_routines.c -lm -O3 -Wall

clean:
	rm -f BiPoS
