
#include <string>
#include <SOFT.h>

/**
 * Prints an error message to stderr.
 * Supports printf syntax.
 *
 * id:  Optional error ID (to only print error once).
 * msg: Message (in printf-format) to print.
 */
template<typename ... Args>
void SOFT::PrintError(const std::string& msg, Args&& ... args) {
    PrintError(SOFT::MESSAGE_GENERAL, msg, std::forward<Args>(args) ...);
}
template<typename ... Args>
void SOFT::PrintError(const SOFT::message_t id, const std::string& msg, Args&& ... args) {
    if (!SOFT::VerifyMessage(id))
        return;

#ifdef COLOR_TERMINAL
    std::string fullmsg = "\x1B[1;31m[ERROR]\x1B[0m "+msg+"\n";
#else
    std::string fullmsg = "[ERROR] "+msg+"\n";
#endif
    fprintf(stderr, fullmsg.c_str(), std::forward<Args>(args) ...);
}

/**
 * Prints a warning message once (if the same warning
 * has already been printed, this function doesn't
 * print anything).
 * Supports printf syntax.
 *
 * id:  Optional warning ID (to only print warning once).
 * msg: Message (in printf-format) to print.
 */
template<typename ... Args>
void SOFT::PrintWarning(const std::string& msg, Args&& ... args) {
    PrintWarning(SOFT::MESSAGE_GENERAL, msg, std::forward<Args>(args) ...);
}
template<typename ... Args>
void SOFT::PrintWarning(const SOFT::message_t id, const std::string& msg, Args&& ... args) {
    if (!SOFT::VerifyMessage(id))
        return;

#ifdef COLOR_TERMINAL
    std::string fullmsg = "\x1B[1;33m[WARNING]\x1B[0m "+msg+"\n";
#else
    std::string fullmsg = "[WARNING] "+msg+"\n";
#endif
    fprintf(stderr, fullmsg.c_str(), std::forward<Args>(args) ...);
}

/**
 * Prints a general information message to stdout.
 * Supports printf syntax.
 *
 * id:  Optional message ID (to only print message once).
 * msg: Message (in printf-format) to print.
 */
template<typename ... Args>
void SOFT::PrintInfo(const std::string& msg, Args&& ... args) {
    PrintInfo(MESSAGE_GENERAL, msg, std::forward<Args>(args) ...);
}
template<typename ... Args>
void SOFT::PrintInfo(const SOFT::message_t id, const std::string& msg, Args&& ... args) {
    if (!SOFT::VerifyMessage(id))
        return;

    std::string fullmsg = msg+"\n";
    fprintf(stdout, fullmsg.c_str(), std::forward<Args>(args) ...);
}

