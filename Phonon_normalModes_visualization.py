import numpy as np
import matplotlib.pyplot as plt
import imageio
import os

def animate_atoms_with_vibrating_arrows(R, modes, atom_info=None, amplitude=0.8, base_arrow_scale=5.5,
                                         frames=60, out_name="vibrating_atoms_with_arrows.gif",
                                         figsize=(10, 8), atom_color='steelblue', arrow_color='crimson'):
    """
    Animate atoms vibrating and arrows oscillating with them, saved as a GIF.
    """
    temp_dir = "_tmp_vib_frames"
    os.makedirs(temp_dir, exist_ok=True)

    for frame in range(frames):
        t = 2 * np.pi * frame / frames
        displacement = amplitude * np.sin(t) * modes
        positions = R + displacement

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_box_aspect([1, 1, 1])
      
        ax.set_axis_off()

        # Set consistent limits with padding
        pad = amplitude * 2
        ax.set_xlim(R[:, 0].min() - pad, R[:, 0].max() + pad)
        ax.set_ylim(R[:, 1].min() - pad, R[:, 1].max() + pad)
        ax.set_zlim(R[:, 2].min() - pad, R[:, 2].max() + pad)

        # Plot atoms
        if atom_info:
            n = len(atom_info[0])
            atom_type = atom_info[0]
            atom_count = atom_info[1]
            pos_color = []
            start = 0
            for c in atom_count:
                end = start + c
                pos_color.append(positions[start:end])
                start = end
            for j in range(n):
                ax.scatter(pos_color[j][:, 0], pos_color[j][:, 1], pos_color[j][:, 2],
                           s=50, label=atom_type[j], edgecolors='k')
            ax.legend()

        else:
            ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                    s=50, color=atom_color, edgecolors='k')

        # Scale arrows with same sine modulation
        arrow_lengths = base_arrow_scale * np.sin(t)
        arrow_vectors = arrow_lengths * modes

        ax.quiver(
            positions[:, 0], positions[:, 1], positions[:, 2],
            arrow_vectors[:, 0], arrow_vectors[:, 1], arrow_vectors[:, 2],
            length=1.0, normalize=False,
            color=arrow_color, arrow_length_ratio=0.2, linewidth=1.5
        )

        # Save frame
        fname = os.path.join(temp_dir, f"frame_{frame:03d}.png")
        plt.savefig(fname, bbox_inches='tight')
        plt.close()

    # Assemble GIF
    images = [imageio.imread(os.path.join(temp_dir, f)) for f in sorted(os.listdir(temp_dir))]
    imageio.mimsave(out_name, images, fps=30, loop=0)

    # Clean up
    for f in os.listdir(temp_dir):
        os.remove(os.path.join(temp_dir, f))
    os.rmdir(temp_dir)

    print(f"✅ Saved animated GIF: {out_name}")

animate_atoms_with_vibrating_arrows(R_gs, modes[881], atom_info = [["B", "N", "C"], [146, 146, 2]], amplitude=0.5, base_arrow_scale=3.4)
