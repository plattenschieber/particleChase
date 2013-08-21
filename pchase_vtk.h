
typedef struct {
        int                 ID;
        double              x[2];
} pchase_particlevtk_t;
typedef struct {
        int                 nParticles;
        pchase_particlevtk_t  *p[25];
}
pchase_quadrant_datavtk_t;
