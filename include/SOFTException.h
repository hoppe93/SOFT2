#ifndef _SOFT_EXCEPTION_H
#define _SOFT_EXCEPTION_H

#include <string>
#include <vector>
#include <softlib/SOFTLibException.h>

/**
 * General-purpose SOFT Exception class.
 */

class SOFTException : public SOFTLibException {
    private:
        std::vector<std::string> modules;
	public:
        SOFTException() {}

		template<typename ... Args>
		SOFTException(const std::string& msg, Args&& ... args) : SOFTLibException(msg, std::forward<Args>(args) ...) {}

        void AddModule(const std::string& m) { modules.push_back(m); }

        std::vector<std::string> &GetModules() { return modules; }
        std::string GetModulesString() {
            std::string m;
            if (modules.size() > 0) {
                std::string m = modules.front();
                for (std::vector<std::string>::iterator it = modules.begin()+1; it != modules.end(); it++) {
                    m += "/" + (*it);
                }

            }

            return m;
        }
};

#endif/*_SOFT_EXCEPTION_H*/
