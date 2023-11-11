#pragma once

/**
 * @brief Face Orientations
 *
 * This enum represents the face orientation flags.
 */
enum class faceorientation {
  p_1, // positive normal direction and edge vertices match with first face edge
  p_2, // positive normal direction and rotated by once
  p_3, // positive normal direction and rotated by twice
  p_4, // positive normal direction and rotated by three times
  n_1, // negative normal direction and edge vertices match with first face edge
  n_2, // negative normal direction and rotated once
  n_3, // negative normal direction and rotated twice
  n_4  // negative normal direction and rotated three times
};

