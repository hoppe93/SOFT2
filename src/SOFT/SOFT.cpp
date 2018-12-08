/**
 * Constructor of SOFT class.
 */

#include "SOFT.h"

bool SOFT::message_checklist[MESSAGE_LAST] = {false};
#ifdef COLOR_TERMINAL
    const std::string SOFT::PRINT_YES="\x1B[1;32mYES\x1B[0m";
    const std::string SOFT::PRINT_NO ="\x1B[1;31mNO\x1B[0m";
#else
    const std::string SOFT::PRINT_YES="YES";
    const std::string SOFT::PRINT_NO ="NO";
#endif

/**
 * Constructor.
 */
SOFT::SOFT() { }
/**
 * Destructor.
 */
SOFT::~SOFT() {}

/**
 * Print a single new-line character in the 'Info' channel.
 */
void SOFT::PrintInfo() {
    PrintInfo(MESSAGE_GENERAL, "");
}

/**
 * Verify against the given checklist that the
 * given message has not been emitted previously.
 *
 * id: Message ID.
 *
 * RETURNS true if the message has NOT been
 * emitted previously (and is safe to emit).
 */
bool SOFT::VerifyMessage(const message_t id) {
    bool retval = false,
         valid = (MESSAGE_GENERAL < id && id < MESSAGE_LAST);

    if (!valid)
        return true;
    else if (message_checklist[id])
        return false;
    
    #pragma omp critical (SOFT_VerifyMessage)
    {
        if (message_checklist[id] == false) {
            message_checklist[id] = true;
            retval = true;
        }
    }

    return retval;
}

