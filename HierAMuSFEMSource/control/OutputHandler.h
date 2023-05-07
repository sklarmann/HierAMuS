// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef OUTPUTHANDLER_H_
#define OUTPUTHANDLER_H_

#include <fstream>
#include <iostream>
#include <memory>
#include <spdlog/common.h>
#include <string>
#include <vector>

#include <iomanip>

#include <filesystem>

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

namespace HierAMuS {

/**
 * @brief Enumerator to handle the output-level from NoLog=0, BasicLog=1 to
 * FullLog=2.
 */
enum class LogLevel {
  NoLog = 0,
  critical = 1,
  error = 2,
  Warning = 3,
  BasicLog = 4,
  Debug = 5,
  FullLog = 6
};

/**
 * @brief Class which handles the output.
 */
class OutputHandler {
public:
  /**
   * @brief Constructor which sets the initial output Level to NoLog.
   */
  OutputHandler();
  ~OutputHandler(){};
  /**
   * @brief Handles the output to the console and the log-file.
   *
   * @tparam T Overload for different objects which support the output stream
   * operator.
   * @param input String to process.
   * @return Returns the OutputHandler.
   */
  template <class T> OutputHandler &operator<<(const T &input) {

    if (this->ConsoleLogLevel >= this->CurrentLogLevel) {
      std::cout << input;
    }

    if (this->FileLogLevel >= this->CurrentLogLevel &&
        this->LogFile.is_open()) {
      this->LogFile << input;
    }

    return *this;
  };
  OutputHandler &operator<<(std::ostream &(*pf)(std::ostream &)) {
    if (this->ConsoleLogLevel >= this->CurrentLogLevel) {
      std::cout << pf;
    }

    if (this->FileLogLevel >= this->CurrentLogLevel &&
        this->LogFile.is_open()) {

      this->LogFile << pf;
    }
    return *this;
  };

  /**
   * @brief Returns the OutputHandler in all state.
   *
   * @return OutputHandler object.
   */
  OutputHandler &all();
  /**
   * @brief Returns the OutputHandler in basic state.
   *
   * @return OutputHandler object.
   */
  OutputHandler &basic();
  /**
   * @brief Returns the OutputHandler in debug state.
   *
   * @return OutputHandler object.
   */
  OutputHandler &debug();

  /**
   * @brief Sets the log-level offset, which gets processed.
   *
   * @param Console Level for console output.
   * @param File Level for log-file output.
   */
  void setLogLevel(LogLevel File, LogLevel Console);
  /**
   * @brief Opens the log-file with name file.
   * @param file String with the filename to open.
   */
  void openLogFile(const std::string &path, const std::string &file);
  /**
   * @brief Closes the log-file.
   */
  void closeLogFile();
  /**
   * @brief Sets the amount of digits to output.
   * @param num Integer containing the amount of digits to output.
   */
  void setPrecision(int num);

  template <typename T>
  void printTwoLineMatrix(std::vector<std::string> names, std::vector<T> values,
                          int fieldWidth, int prec) {
    int tempPrec = this->numberWidth;
    this->setPrecision(prec);
    *this << "   ";
    for (auto i = names.begin(); i != names.end(); ++i) {
      *this << std::left << std::setw(fieldWidth) << std::setfill(' ') << *i;
    }
    *this << std::endl << "   ";
    for (auto i = values.begin(); i != values.end(); ++i) {
      *this << std::left << std::setw(fieldWidth) << std::setfill(' ') << *i;
    }
    *this << std::endl;

    this->setPrecision(tempPrec);
  }

  int getNumberWidth() { return this->numberWidth + 9; };
  int getNumberPrecision() { return this->numberWidth; };

  void output(LogLevel print, LogLevel write) {
    if (print <= this->ConsoleLogLevel) {
      std::cout << std::endl;
    }
    if (write <= this->FileLogLevel) {
      if (this->LogFile.is_open()) {
        this->LogFile << std::endl;
      }
    }
  }

  template <typename T, typename... Args>
  void output(LogLevel print, LogLevel write, T firstArg,
              Args... remainingArgs) {
    if (print <= this->ConsoleLogLevel) {
      std::cout << firstArg;
    }
    if (write <= this->FileLogLevel) {
      if (this->LogFile.is_open()) {
        this->LogFile << firstArg;
      }
    }
    this->output(print, write, remainingArgs...);
  }

  void outputString(LogLevel print, LogLevel write, std::string toPrint) {

    auto ll = print < write ? print : write;

    if (ll == LogLevel::NoLog) {

    } else if (ll == LogLevel::critical) {
      this->Logger->critical(toPrint);
    } else if (ll == LogLevel::error) {
      this->Logger->error(toPrint);
    } else if (ll == LogLevel::Warning) {
      this->Logger->warn(toPrint);
    } else if (ll == LogLevel::BasicLog) {
      this->Logger->info(toPrint);
    } else if (ll == LogLevel::Debug) {
      this->Logger->debug(toPrint);
    } else if (ll == LogLevel::FullLog) {
      this->Logger->trace(toPrint);
    }
  }

  void outputStringBasic(std::string toPrint) { this->Logger->info(toPrint); }
  void outputStringDebug(std::string toPrint) { this->Logger->debug(toPrint); }

  auto getSPDLogger() -> spdlog::logger & { return *this->Logger; }

protected:
  LogLevel CurrentLogLevel;               // Outputloglevel
  LogLevel FileLogLevel, ConsoleLogLevel; // Loglevel

  std::ofstream LogFile;
  int numberWidth = 6;

private:
  std::shared_ptr<spdlog::sinks::basic_file_sink_mt> fileSink;
  std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> consoleSink;
  std::shared_ptr<spdlog::logger> Logger;
};

} /* namespace HierAMuS */

#endif /* OUTPUTHANDLER_H_ */
