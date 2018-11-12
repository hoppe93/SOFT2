/**
 * Implementation of the ToolHandler class.
 *
 * This class keeps a list of each enabled tool (with one handler
 * per thread) and allows the SOFTLocal object to invoke all tools
 * with only one call to the ToolHandler object.
 */

#include <softlib/config.h>
#include "Orbit/Orbit.h"
#include "PhaseSpace/Particle.h"
#include "Tools/Tool.h"
#include "Tools/ToolHandler.h"

/**
 * Constructor
 */
ToolHandler::ToolHandler() { }

/**
 * Adds the given tool to the list of tools
 * of this handler. Note that a direct reference
 * to the Tool is added, and thus it is _not_
 * just copied.
 * 
 * t: Pointer to Tool to add.
 */
void ToolHandler::AddTool(Tool *t) {
    this->tools.push_back(t);
}

/**
 * Calls the 'Finish' method on each of the tools
 * handled by this ToolHandler object.
 */
void ToolHandler::Finish() {
    for (size_t i = 0; i < tools.size(); i++) {
        tools[i]->Finish();
    }
}

/**
 * Calls the 'Handle' method on each of the tools
 * handled by this ToolHandler object, passing on
 * the given Orbit object.
 *
 * o: Orbit object to handle.
 */
void ToolHandler::Handle(Orbit *o, Particle *p) {
    for (size_t i = 0; i < tools.size(); i++) {
        tools[i]->Handle(o, p);
    }
}

/**
 * Calls the 'Output' method on each of the tools
 * handled by this ToolHandler object.
 */
void ToolHandler::Output() {
    for (size_t i = 0; i < tools.size(); i++) {
        tools[i]->Output();
    }
}

