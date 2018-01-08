#version 120

varying vec2 v_sky_texture_coord;

uniform sampler2D u_sky_texture;

void main() {
    vec3 rgb = texture2D(u_sky_texture, v_sky_texture_coord).rgb;

    gl_FragColor.rgb = clamp(rgb,0.0,1.0);
    gl_FragColor.a = 1;
}