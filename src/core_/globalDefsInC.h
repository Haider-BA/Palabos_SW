



#ifndef GLOBALDEFSINC_H_INCLUDED
#define GLOBALDEFSINC_H_INCLUDED


struct UniformSequence
{
    int num;
    int length;
    void* sequence;
};

struct UniformSequence_const
{
    int num;
    int length;
    const void* sequence;
};


struct NonUniformSequence
{
    int num;
    int* bar;
    void* sequence;
};

struct NonUniformSequence_const
{
    int num;
    int* bar;
    const void* sequence;
};

struct KeySequence
{
    int num;
    int* length;
    int* position;
    void *sequence;
};

struct KeySequence_const
{
    int num;
    int* length;
    int* position;
    const void *sequence;
};


#endif // GLOBALDEFSINC_H_INCLUDED
