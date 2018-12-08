#ifndef _TOOL_HANDLER_H
#define _TOOL_HANDLER_H

#include <vector>
#include <softlib/config.h>
#include "Orbit/Orbit.h"
#include "PhaseSpace/Particle.h"
#include "Tools/Tool.h"

class ToolHandler {
    private:
        std::vector<Tool*> tools;
    public:
        ToolHandler();

        void AddTool(Tool*);
        void Finish();
        void Handle(Orbit*, Particle*);
        void Initialize();
        void Output();

        const std::vector<Tool*> &GetTools() const { return tools; }
};

#endif/*_TOOL_HANDLER_H*/
