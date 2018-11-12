#ifndef _TOOL_H
#define _TOOL_H

#include <string>
#include <softlib/Configuration.h>
#include "Orbit/Orbit.h"
#include "PhaseSpace/Particle.h"

class Tool {
	private:
		std::string name;
        std::string type;
	
	public:
		Tool(const std::string);

        virtual void Configure(struct global_settings*, ConfigBlock*, ConfigBlock*) = 0;
        virtual void Finish() = 0;
        virtual void Handle(Orbit*, Particle*) = 0;
        virtual void Output() = 0;
        virtual void Welcome(const std::string &prefix="  ") = 0;

        const std::string GetName() const { return name; }
        const std::string GetType() const { return type; }
        const std::string GetTypeAndName() const { return ("@"+type+"("+name+")"); }
        void SetName(const std::string &name) { this->name = name; }
};

class ToolException : public SOFTException {
    public:
        template<typename ... Args>
        ToolException(const std::string &msg, Args&& ... args)
            : SOFTException(msg, std::forward<Args>(args) ...) {
            AddModule("Tool");
        }
};

#endif/*_TOOL_H*/
