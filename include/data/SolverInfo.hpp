#ifndef SOLVERINFO_HPP
#define SOLVERINFO_HPP

#include "Traits.hpp"
#include <map>

class SolverInfo {
public:
    /**
     * Typedef for solver result
     */
    typedef std::map<std::string, Real> result_t;

    enum TableLines { NO_TABLE_LINES = -1 };


    const bool exportIterationVTK;
    const int stepOutputMaxNumber;
    const std::string VTKfinalSolutionPath;
    const std::string DATfinalSolutionPath;
    const Vector DATextractionPoint;
    result_t results;

    SolverInfo(const bool exportIterationVTK,
               const int stepOutputMaxNumber,
               const std::string &VTKfinalSolutionPath,
               const std::string &DATfinalSolutionPath,
               const Vector &DATextractionPoint) : exportIterationVTK(exportIterationVTK),
                                                   stepOutputMaxNumber(stepOutputMaxNumber > 0 ? stepOutputMaxNumber : -1),
                                                   VTKfinalSolutionPath(VTKfinalSolutionPath),
                                                   DATfinalSolutionPath(DATfinalSolutionPath),
                                                   DATextractionPoint(DATextractionPoint) {
    };
};

#endif //SOLVERINFO_HPP
