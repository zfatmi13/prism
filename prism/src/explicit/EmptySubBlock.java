/*
 * Copyright (C) 2020  Zainab Fatmi and Franck van Breugel
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You can find a copy of the GNU General Public License at
 * <http://www.gnu.org/licenses/>.
 */

package explicit;

/**
 * An empty sub-block.
 * 
 * @author Zainab Fatmi
 * @author Franck van Breugel
 */
public final class EmptySubBlock extends SubBlock {

	/**
	 * Tests whether the given object is an empty sub-block.
	 * 
	 * @param object an object
	 * @return true if the given object is an empty sub-block, false otherwise
	 */
	@Override
	public boolean equals(Object object) {
		return object instanceof EmptySubBlock;
	}

	/**
	 * Throws an exception since the states in this empty sub-block do not
	 * transition to any of the blocks that have been split in the previous
	 * refinement step.
	 * 
	 * @throws UnsupportedOperationException always
	 */
	@Override
	public Double get(Integer block) {
		throw new UnsupportedOperationException();
	}

	/**
	 * Throws an exception since the states in this empty sub-block do not
	 * transition to any of the blocks that have been split in the previous
	 * refinement step.
	 * 
	 * @throws UnsupportedOperationException always
	 */
	@Override
	public void put(Integer block, Double probability) {
		throw new UnsupportedOperationException();
	}
}
