// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <iostream>
#include <memory>
#include <fstream>
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
  void openLogFile(const std::string &path, const std::string &file, bool override = true);
  /**
   * @brief Closes the log-file.
   */
  void closeLogFile();
  /**
   * @brief Sets the amount of digits to output.
   * @param num Integer containing the amount of digits to output.
   */
  void setPrecision(int num);




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

private:
  std::shared_ptr<spdlog::sinks::basic_file_sink_mt> fileSink;
  std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> consoleSink;
  std::shared_ptr<spdlog::logger> Logger;
};

} /* namespace HierAMuS */

