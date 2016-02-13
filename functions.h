#include <gsl/gsl_rng.h>

class Generator{
	public:
		Generator();
		Generator(size_t s);
		~Generator();
		void set_seed(size_t seed);
		size_t seed;
		size_t count;
		const gsl_rng_type * T;
		gsl_rng * r;
		double uniform();

};

class Fiber{
	public:
		Fiber();
		Fiber(double x, double y, double z, double d1, double d2, double n_x, double n_y, double n_z, unsigned int N);
		void generateSources(Generator & generator);

	private:
		double x; //center of fiber
		double y;
		double z;

		double d1; //core diameter
		double d2; //cladding diameter

		double n_x; //Surface normal directional cosines
		double n_y;
		double n_z;

		unsigned int N; //Number of source points
};


