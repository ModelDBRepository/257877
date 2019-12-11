This is the model associated with the paper:

Assisi Stopfer Bazhenov

The code was contributed by Collins Assisi. This model is written in
C++. The following instructions work have worked under Ubuntu 13.04
although the code should work with little or no modification on other
platforms with ansi compilers.

You can use the supplied script to compile and run which you can
execute on the linux command line

. submitJobAlphaSpikingSTDP2.sh
or you can run each step by hand:

How to compile:

Change to the directory in which you have unziped the archive file,
then type:

g++ KCAlphaSpikingSTDP2.cpp -lm -O4 -ffast-math -o n

How to run:
./n

The code models the response of Mushroom Body output neurons receiving
convergent input from Kenyon cells. MBONs were modeled as

Rulkov-Bazhenov abstract map based neurons (Rulkov NF, Timofeev I,
Bazhenov M (2004) Oscillations in large-scale cortical networks:
map-based model. J Comput Neurosci 17:203-23 ).

The synaptic strength evolved in an activity dependent
manner. Synapses were modeled according to

Jesper Sjöström and Wulfram Gerstner (2010), Scholarpedia, 5(2):1362.

Training data is given in spikesKC2. This data is read by
KCAlphaSpikingSTDP2.cpp and used as input to set of MBONs.
