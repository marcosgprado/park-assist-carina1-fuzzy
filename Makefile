all:

	#gcc teste1.c -o teste1 -I/usr/local/include/player-3.0 -L/usr/local/lib -lplayerc
	#gcc teste2.c -o teste2 -I/usr/local/include/player-3.0 -L/usr/local/lib -lplayerc
	g++ -c main.cpp -I/usr/local/include/player-3.0
	g++ -c main01.cpp -I/usr/local/include/player-3.0
	g++ -c main2_0.cpp -I/usr/local/include/player-3.0
	g++ -c main3_0.cpp -I/usr/local/include/player-3.0
	g++ -c fis.cpp -I/usr/local/include/player-3.0

	g++ main.o fis.o -o main -L/usr/local/lib `pkg-config --cflags --libs playerc++` -lpthread -pipe  -W -Wall -O2
	g++ main01.o -o main01 -L/usr/local/lib `pkg-config --cflags --libs playerc++` -lpthread -pipe  -W -Wall -O2
	g++ main2_0.o fis.o -o main2_0 -L/usr/local/lib `pkg-config --cflags --libs playerc++` -lpthread -pipe  -W -Wall -O2
	g++ main3_0.o fis.o -o main3_0 -L/usr/local/lib `pkg-config --cflags --libs playerc++` -lpthread -pipe  -W -Wall -O2

	
