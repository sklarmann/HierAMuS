// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <iostream>
#include <fstream>
#include "OutputHandler.h"

namespace HierAMuS {
  
  OutputHandler::OutputHandler() {
    CurrentLogLevel = LogLevel::NoLog;
    FileLogLevel = LogLevel::NoLog;
    ConsoleLogLevel = LogLevel::NoLog;
  };

  auto OutputHandler::all() -> OutputHandler& {
    this->CurrentLogLevel = LogLevel::NoLog;
    return *this;
  }

  
  auto OutputHandler::basic() -> OutputHandler& {
    this->CurrentLogLevel = LogLevel::BasicLog;
    return *this;
  }

  
  auto OutputHandler::debug() -> OutputHandler& {
    this->CurrentLogLevel = LogLevel::FullLog;
    return *this;
  }

  void OutputHandler::setLogLevel(LogLevel File, LogLevel Console){
    this->ConsoleLogLevel = Console;
    this->FileLogLevel = File;

    auto switcher = [](LogLevel lvl){
      switch (lvl) {
        case LogLevel::NoLog:
          return spdlog::level::off;
        case LogLevel::critical:
          return spdlog::level::critical;
          break;
        case LogLevel::error:
          return spdlog::level::err;
          break;
        case LogLevel::Warning:
          return spdlog::level::warn;
          break;
        case LogLevel::BasicLog:
          return spdlog::level::info;
          break;
        case LogLevel::Debug:
          return spdlog::level::debug;
          break;
        case LogLevel::FullLog:
          return spdlog::level::trace;
          break;
        default:
          return spdlog::level::info;
      
      }
    };

    this->consoleSink->set_level(switcher(Console));
    this->fileSink->set_level(switcher(File));


  }


  void OutputHandler::openLogFile(const std::string &path, const std::string &file, bool override){

    std::filesystem::path directory(path);
    std::filesystem::path logfilename(file + ".log");


    
    this->consoleSink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    this->consoleSink->set_level(spdlog::level::info);

    std::filesystem::path spdLogFile =
        directory / std::filesystem::path(file + ".spdlog");
    this->fileSink =
        std::make_shared<spdlog::sinks::basic_file_sink_mt>(spdLogFile.string(),override);
    this->fileSink->set_level(spdlog::level::info);

    std::vector<spdlog::sink_ptr> sinks{consoleSink, fileSink};

    this->Logger = std::make_shared<spdlog::logger>(file, sinks.begin(), sinks.end());
    this->Logger->set_level(spdlog::level::trace);
  }

  void OutputHandler::closeLogFile() {

  }

} /* namespace HierAMuS */
