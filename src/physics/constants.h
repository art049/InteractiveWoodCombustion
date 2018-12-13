extern const float Deltat[1] {0.0001f}; // Easier to pass array to cuda
extern const uint GRID_WIDTH =  480;
extern const uint GRID_HEIGHT=  480;
extern const uint GRID_DEPTH =  288;
extern const dim3 M_i { 16 , 16 , 4  };