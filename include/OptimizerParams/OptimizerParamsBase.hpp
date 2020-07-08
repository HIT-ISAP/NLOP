#ifndef OPTIMIZERPARAMSBASE_HPP
#define OPTIMIZERPARAMSBASE_HPP

#include <Utils/Functor.hpp>
#include <string>

namespace NLOP {
/// @class NLOP::OptimizerParamsBase
/// @brief Abstract base class for params of all non-linear optimization methods
class OptimizerParamsBase
{
public:
    // How much information of the optimization process expected to be printed
    enum Verbosity{
        SILENT,  // No information printed
        SUMMARY, // Information of the initialization and termination printed
        DETAIL   // Information of every iteration step in the optimization process printed
    };

    /// @brief Verbosity translator: from Verbosity to std::string
    std::string verbosityTranslator(Verbosity value)
    {
        std::string s;
        switch (value) {
        case SILENT:
            s = "SILENT";
            break;
        case SUMMARY:
            s = "SUMMARY";
            break;
        case DETAIL:
            s = "DETAIL";
            break;
        default:
            s = "UNDEFINED";
            break;
        }
        return s;
    }

    /// @brief Verbosity translator: from std::string to Verbosity
    Verbosity verbosityTranslator(const std::string& value)
    {
        std::string s = value;
        if (s == "SILENT")
            return SILENT;
        if (s == "SUMMARY")
            return SUMMARY;
        if (s == "DETAIL")
            return DETAIL;

        return SUMMARY; // If invalid value is received, choose SUMMARY mode
    }

    /// @brief Use default optimizer params
    virtual void setDefaults()
    {
        iteration_times = 0;
        max_iteration_times = 10000;
        verbosity = SUMMARY;
        log_file = false;
    }

    /// @brief Get and set max iteration times
    void setMaxItertaions(size_t value) { max_iteration_times = value; }
    size_t getMaxIterations() const { return max_iteration_times; }

    /// @brief Get and set verbosity
    void setVerbosity(std::string value) { verbosity = verbosityTranslator(value); }
    Verbosity getVerbosity() { return verbosity; }

    /// @brief Get and set recent iteration times
    void setIterationTimes(size_t value) { iteration_times = value; }
    size_t getIterationTimes() const { return iteration_times; }

    /// @brief Enable and disable logging file
    void enableLog() { log_file = true; }
    void disableLog() { log_file = false; }

    bool isLogFile() const { return log_file; }

    /// @brief Recent iteration times + 1
    void nextIteration() { iteration_times++; }

    /// @brief Print optimizer params
    virtual void print()
    {
        std::cout << "maximum iterations: " << max_iteration_times << "\n";;
        std::cout << "verbosity: " << verbosityTranslator(verbosity) << std::endl;
    }

protected:
    size_t iteration_times;     // Record iteration times in the optimization process (initial = 0)
    size_t max_iteration_times; // Maximum iteration times to stop iterating (default 10000)
    Verbosity verbosity;        // How much information of the optimization process expected to be printed (default = SUMMARY)
    bool log_file;              // Whether to write log file

    virtual ~OptimizerParamsBase() {}
};
}


#endif
