/**
 *
 * @file   femfussion.h
 * @brief  Main file of the FEMFFUSION program.
 *
 */

#ifndef FEMFUSSION_H
#define FEMFUSSION_H

#include <deal.II/base/parameter_handler.h>
#include "utils.h"

/**
 * @brief Print the FEMFUSION logo.
 */
void print_logo(std::ostream &out);

/**
 *  @brief Declare all the entries defined in the prm file.
 */
void prm_declare_entries(ParameterHandler &prm);


/**
 *  @brief Declare the entries defined in the dynamic prm file.
 */
void prm_dyn_entries(ParameterHandler &prm);


/**
 * Main function.
 */
int main(int argc, char **argv);

#endif /* FEMFUSSION_H */
