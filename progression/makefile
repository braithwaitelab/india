model: hiv_toxicity.o distributions.o perlin.o constantsHeader.o
	g++ -O3 -o model hiv_toxicity.o distributions.o perlin.o constantsHeader.o

hiv_toxicity.o:  hiv_toxicity.cc hiv.h distributions.h perlin.h constantsHeader.h output.h
	g++ -O3 -c hiv_toxicity.cc

distributions.o:  distributions.cc distributions.h
	g++ -O3 -c distributions.cc

perlin.o:  perlin.cc perlin.h
	g++ -O3 -c perlin.cc

constantsHeader.o:  constantsHeader.cc constantsHeader.h
	g++ -O3 -c constantsHeader.cc

clean:
	rm model hiv_toxicity.o distributions.o perlin.o constantsHeader.o
