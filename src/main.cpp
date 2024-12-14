#include <spdlog/cfg/env.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <argparse/argparse.hpp>

#include "version.h"
#include "markdup/md_args.h"

using namespace std;

extern MarkDupsArg gMdArg;

int main(int argc, char *argv[]) {
    // init log
    spdlog::set_default_logger(spdlog::stderr_color_st("fastdup"));
    spdlog::cfg::load_env_levels();
    // init arg parser
    argparse::ArgumentParser cliArgs()


    spdlog::info("fast markduplicates start");

    return 0;
}