#include "utils/logger.h"
#include "gds/gas/types.h"
#include "gds/gas/corrector.h"

int main(int argc, char * argv[])
{
    Logger logger;
    LOG("GDS online");
    gds::SimpleCell<double> l;
    gds::SimpleCell<double> m;
    gds::SimpleCell<double> r;
    gds::corrector(l, m, r, 0.1);
    return 0;
}
