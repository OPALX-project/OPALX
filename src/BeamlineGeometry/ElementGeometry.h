#ifndef OPALX_ElementGeometry_HH
#define OPALX_ElementGeometry_HH

#include "BeamlineGeometry/Port.h"

/**
 * @class ElementGeometry
 * @brief Named body-relative ports used for placement and chaining.
 *
 * This is a lightweight implementation holder for the canonical local chart and
 * its named interface frames. In the language of the placement note, these
 * ports represent marked boundary charts with adapted frames. The first
 * redesign stage keeps the set minimal: `entry`, `body`, and `exit`.
 */
class ElementGeometry {
public:
    ElementGeometry()
        : entry_m("entry", CoordinateSystemTrafo()),
          body_m("body", CoordinateSystemTrafo()),
          exit_m("exit", CoordinateSystemTrafo()) {}

    ElementGeometry(const Port& entry, const Port& body, const Port& exit)
        : entry_m(entry), body_m(body), exit_m(exit) {}

    const Port& getEntry() const { return entry_m; }
    const Port& getBody() const { return body_m; }
    const Port& getExit() const { return exit_m; }

private:
    Port entry_m;
    Port body_m;
    Port exit_m;
};

#endif
