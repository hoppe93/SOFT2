#ifndef _SMPI_H
#define _SMPI_H

#ifndef WITH_MPI
#   error "The SMPI module requires SOFT to be compiled with MPI."
#endif

#include <mpi.h>
#include "SOFTException.h"

namespace SMPI {
    extern MPI_Datatype MPI_SLIBREAL_T;

    void init(int*, char***);
    void finalize();
    void verify(int);

    class MPIException : public SOFTException {
        private:
            int error_code;
        public:
            template<typename ... Args>
            MPIException(const std::string &msg, Args&& ... args)
                : SOFTException(msg, std::forward<Args>(args) ...) {
                AddModule("SMPI");
            }
            MPIException(int err) {
                AddModule("SMPI");
                this->error_code = err;

                char str[MPI_MAX_ERROR_STRING];
                int len;
                MPI_Error_string(err, str, &len);

                ConstructErrorMessage("MPI Error %d: %s.", err, str);
            }

            int GetErrorCode() { return this->error_code; }
    };

}

#endif/*_SMPI_H*/
