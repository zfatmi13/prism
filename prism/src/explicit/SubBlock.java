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
 * A sub-block of a block. During the refinement step, blocks are split into
 * sub-blocks. A sub-block has an ID, which will be used in the next refinement
 * step when the sub-block plays the role of a block.
 * 
 * @author Zainab Fatmi
 * @author Franck van Breugel
 */
public abstract class SubBlock {

	/**
	 * The ID of this sub-block. It will be used in the next refinement step when
	 * the sub-block plays the role of a block
	 */
	private int ID;

	/**
	 * Returns the ID of this sub-block.
	 * 
	 * @return the ID of this sub-block
	 */
	public int getID() {
		return this.ID;
	}

	/**
	 * Sets the ID of this sub-block to the given ID.
	 * 
	 * @param ID the new ID of this sub-block
	 */
	public void setID(int ID) {
		this.ID = ID;
	}

	/**
	 * Returns the probability of transitioning to states in the given block.
	 * 
	 * @param block the ID of a block
	 * @return the probability of transitioning to states in the given block
	 */
	public abstract Double get(Integer block);

	/**
	 * Sets the probability of transitioning to states in the given block to the
	 * given probability.
	 * 
	 * @param block       the ID of a block
	 * @param probability the probability of transitioning to states in the block
	 */
	public abstract void put(Integer block, Double probability);
}
