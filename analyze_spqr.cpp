#include <stdio.h>
#include <fstream>

int main(int argc, char * argv[])
{
    long n, fm, fn, fp;

    /// Getting the number of rows, columns
    std::ifstream f(argv[1]);
    while(!f.eof())
    {
        f >> n >> fm >> fn >> fp;
        if (fm > 10000 || fp > 100)
           printf("%ld: fm %ld fn %ld fp %ld\n", n, fm, fn, fp);
    }
    f.close();
}
