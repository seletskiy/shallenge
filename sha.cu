// compile:
//  v=$(nvidia-smi --query-gpu=compute_cap --format=csv | tail -n1 | tr -cd
//  '[[:digit:]]') nvcc -arch compute_$v -code sm_$v sha.cu -o sha
#include <stdint.h>
#include <stdio.h>
#include <time.h>

#define SHR(a, b) (__funnelshift_r((a), 0, (b)))
#define ROT(a, b) (__funnelshift_r((a), (a), (b)))
#define CH(x, y, z) (((x) & ((y) ^ (z))) ^ (z))
#define MAJ(x, y, z) (((x) & ((y) | (z))) | ((y) & (z)))
#define EP0(x) (ROT((x), 2) ^ ROT((x), 13) ^ ROT((x), 22))
#define EP1(x) (ROT((x), 6) ^ ROT((x), 11) ^ ROT((x), 25))
#define SIG0(x) (ROT((x), 7) ^ ROT((x), 18) ^ SHR((x), 3))
#define SIG1(x) (ROT((x), 17) ^ ROT((x), 19) ^ SHR((x), 10))

#define K                                                                      \
  {0x428A2F98, 0x71374491, 0xB5C0FBCF, 0xE9B5DBA5, 0x3956C25B, 0x59F111F1,     \
   0x923F82A4, 0xAB1C5ED5, 0xD807AA98, 0x12835B01, 0x243185BE, 0x550C7DC3,     \
   0x72BE5D74, 0x80DEB1FE, 0x9BDC06A7, 0xC19BF174, 0xE49B69C1, 0xEFBE4786,     \
   0x0FC19DC6, 0x240CA1CC, 0x2DE92C6F, 0x4A7484AA, 0x5CB0A9DC, 0x76F988DA,     \
   0x983E5152, 0xA831C66D, 0xB00327C8, 0xBF597FC7, 0xC6E00BF3, 0xD5A79147,     \
   0x06CA6351, 0x14292967, 0x27B70A85, 0x2E1B2138, 0x4D2C6DFC, 0x53380D13,     \
   0x650A7354, 0x766A0ABB, 0x81C2C92E, 0x92722C85, 0xA2BFE8A1, 0xA81A664B,     \
   0xC24B8B70, 0xC76C51A3, 0xD192E819, 0xD6990624, 0xF40E3585, 0x106AA070,     \
   0x19A4C116, 0x1E376C08, 0x2748774C, 0x34B0BCB5, 0x391C0CB3, 0x4ED8AA4A,     \
   0x5B9CCA4F, 0x682E6FF3, 0x748F82EE, 0x78A5636F, 0x84C87814, 0x8CC70208,     \
   0x90BEFFFA, 0xA4506CEB, 0xBEF9A3F7, 0xC67178F2}

#define B64(x) ((x) + ((x) > 11 ? 53 + ((x) > 37) * 6 : (46 - (!(x)) * 3)))

#define W(a, b, c, d) d, c, b, a
#define MSG                                                                    \
  W('s', 'e', 'l', 'e'), W('t', 's', 'k', 'i'), W('y', '/', 'H', 'i'),         \
      W('r', 'e', 'M', 'e'), W('/', 'H', 'i', 'H'), W('N', '/', 'H', 'i'),     \
      W('A', 'l', 'e', 'x'), W('/', '2', '5', 'G'), W('H', 's', '/', 'R'),     \
      W('T', 'X', '4', '0'), W('9', '0', '/', '0'), 0, '_', '_', '_', 0, '_',  \
      '_', '_', 0x80, 0, 0, 0, 0x00, 0x00, 0x00, 0x00, 0xb8, 0x01, 0x00, 0x00, \
                                                                               \
      0x3c, 0x77, 0x7a, 0xce, 0x6f, 0x84, 0x5c, 0xbe

__global__ void kernel(uint32_t _a, uint32_t _b, uint32_t _c, uint32_t _d,
                       uint32_t _e, uint32_t _f, uint32_t _g, uint32_t _h,
                       uint32_t _m_18, uint32_t batch, uint64_t chunk,
                       uint8_t *block, uint32_t *target, uint32_t *mutex) {
  uint8_t msg[32 * 16] = {MSG};

  static const __constant__ uint32_t k[64] = K;

  *((uint32_t *)&msg[44]) = batch;
  msg[48] = blockIdx.x;
  msg[49] = blockIdx.y;
  msg[50] = blockIdx.z;
  msg[51] = threadIdx.x;
  msg[53] = threadIdx.y;
  *((uint64_t *)&msg[48]) += chunk;

  uint64_t target_a = ~0;

  for (uint8_t u = 0x41; u < 0x59; u++) {
    for (uint8_t v = 0x41; v < 0x59; v++) {
      msg[54] = u;
      msg[55] = v;

      uint32_t a, b, c, d, e, f, g, h, t1, t2;

      a = _a;
      b = _b;
      c = _c;
      d = _d;
      e = _e;
      f = _f;
      g = _g;
      h = _h;

      uint32_t i = 0, *m = (uint32_t *)msg;

#pragma unroll 4
      for (i = 12; i < 16; ++i) {
        t1 = h + EP1(e) + CH(e, f, g) + k[i] + m[i];
        t2 = EP0(a) + MAJ(a, b, c);
        h = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
        //if (i == 10) { printf("a = 0x%08x;\n", a); printf("b = 0x%08x;\n", b); printf("c = 0x%08x;\n", c); \
        printf("d = 0x%08x;\n", d); printf("e = 0x%08x;\n", e); printf("f = 0x%08x;\n", f); \
        printf("g = 0x%08x;\n", g); printf("h = 0x%08x;\n", h); }
      }

      /* m[16] = _m_16; */
      /* m[17] = _m_17; */
      /* m[16] = 0xce7a773c; */
      /* m[17] = 0xbe5c846f; */
      m[18] = _m_18;

#pragma unroll 43
      for (i = 19; i < 62; ++i) {
        m[i] = SIG1(m[i - 2]) + m[i - 7] + SIG0(m[i - 15]) + m[i - 16];
        /* if (i==18) printf("%08x\n", m[i]); */
        // if (i < 18) { printf("m[%d] = 0x%08x;\n", i, m[i]); }
      }

#pragma unroll 44
      for (i = 16; i < 60; ++i) {
        t1 = h + EP1(e) + CH(e, f, g) + k[i] + m[i];
        t2 = EP0(a) + MAJ(a, b, c);
        h = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
      }

      // unroll four iterations to avoid writing to m[i] and other extra
      // variables
      t1 = h + EP1(e) + CH(e, f, g) + k[60] + m[60];

      h = d + t1;
      d = t1 + EP0(a) + MAJ(a, b, c);

      g = g + EP1(h) + CH(h, e, f) + k[61] + m[61];
      c = c + g;

      t1 = f + EP1(c) + CH(c, h, e) + k[62] + SIG1(m[60]) + m[55] +
           SIG0(m[47]) + m[46];
      t2 = EP0(d) + MAJ(d, a, b) + g;

      f = t1 + EP0(t2) + MAJ(t2, d, a);
      // we are interested in checking second part of the hash first
      if (uint64_t(f + 0xbb67ae85) < *target) {
        t1 = t1 + b;
        b = e + EP1(t1) + CH(t1, c, h) + k[63] + SIG1(m[61]) + m[56] +
            SIG0(m[48]) + m[47] + EP0(f) + MAJ(f, t2, d);

        if (!uint64_t(b - 0x95f61999)) {
          target_a =
              (uint64_t(f + 0xbb67ae85)) + (uint64_t(b - 0x95f61999) << 32);
          goto found;
        }
      }
      // if (target_a < *target) goto found;
    }
  }
found:

  if (target_a > *target)
    return;

  while (atomicCAS(mutex, 0, 1) != 0)
    ;
  if (target_a < *target) {
    *target = target_a;
    block[3] = msg[44];
    block[2] = msg[45];
    block[1] = msg[46];
    block[0] = msg[47];

    block[7] = msg[48];
    block[6] = msg[49];
    block[5] = msg[50];
    block[4] = msg[51];

    block[11] = 0;
    block[10] = msg[53];
    block[9] = msg[54];
    block[8] = msg[55];
  }
  atomicExch(mutex, 0);
}

#define cudaNoError(expr)                                                      \
  { cudaAssert((expr), __FILE__, __LINE__); }
inline void cudaAssert(cudaError_t code, const char *file, int line) {
  if (code != cudaSuccess) {
    fprintf(stderr, "cudaAssert: %s %s:%d\n", cudaGetErrorString(code), file,
            line);
    exit(code);
  }
}

int main(int argc, char *argv[]) {
  setvbuf(stdout, NULL, _IONBF, 0);

  uint32_t shard;
  uint32_t target, batch;

  if (argc < 3)
    fprintf(stderr, "usage: %s <shard> <target> <batch start>\n", argv[0]),
        exit(1);
  if (sscanf(argv[1], "%x", &shard) == EOF)
    exit(2);
  if (sscanf(argv[2], "%x", &target) == EOF)
    exit(2);
  if (sscanf(argv[3], "%x", &batch) == EOF)
    exit(2);

  fprintf(stderr, "* shard %02x | target start %08x | batch start %06x\n",
          shard, target, batch);

  uint8_t block[12] = {0}, *block_cu;

  const uint8_t chunks_n = 3;
  uint32_t *target_cu[chunks_n], *mutex_cu[chunks_n];

  cudaNoError(cudaMalloc(&block_cu, sizeof(uint8_t) * 12));

  for (uint8_t j = 0; j < chunks_n; j++) {
    cudaNoError(cudaMalloc(&target_cu[j], sizeof(uint32_t)));
    cudaNoError(cudaMalloc(&mutex_cu[j], sizeof(uint32_t)));
  }

  for (uint8_t j = 0; j < chunks_n; j++) {
    cudaNoError(cudaMemcpy(target_cu[j], &target, sizeof(uint32_t),
                           cudaMemcpyHostToDevice));
    cudaNoError(cudaMemset(mutex_cu[j], 0, sizeof(uint32_t)));
  }

  cudaNoError(cudaMemset(block_cu, 0, sizeof(uint8_t) * 12));

  const uint64_t batches_n = 64 * 64 * 64;

  static const uint32_t k[64] = K;

  uint32_t a_0 = 0x6a09e667;
  uint32_t b_0 = 0xbb67ae85;
  uint32_t c_0 = 0x3c6ef372;
  uint32_t d_0 = 0xa54ff53a;
  uint32_t e_0 = 0x510e527f;
  uint32_t f_0 = 0x9b05688c;
  uint32_t g_0 = 0x1f83d9ab;
  uint32_t h_0 = 0x5be0cd19;

  /* #undef W */
  /* #define W(a, b, c, d) a, b, c, d */
  uint8_t msg[64 + 4 * 4] = {MSG};
  /* #undef W */

  for (; batch < batches_n; batch++) {
    struct timespec t_start, t_end;
    uint64_t batch_size = 0;

    msg[44] = B64((batch >> 0) & 0x3f);
    msg[45] = B64((batch >> 6) & 0x3f);
    msg[46] = B64((batch >> 12) & 0x3f);
    msg[47] = B64(shard);

#undef ROT
#undef SHR
#define ROT(a, b) (((a) >> (b)) | ((a) << (32 - (b))))
#define SHR(a, b) ((a) >> (b))

    uint32_t a = a_0;
    uint32_t b = b_0;
    uint32_t c = c_0;
    uint32_t d = d_0;
    uint32_t e = e_0;
    uint32_t f = f_0;
    uint32_t g = g_0;
    uint32_t h = h_0;
    for (uint8_t i = 0; i < 12; ++i) {
      uint32_t t1 = h + EP1(e) + CH(e, f, g) + k[i] + ((uint32_t *)(msg))[i];
      uint32_t t2 = EP0(a) + MAJ(a, b, c);
      h = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
    }

    uint32_t *m = (uint32_t *)msg;
    /* m[16] = 0xce7a773c; */
    /* m[17] = 0xbe5c846f; */
    for (uint32_t i = 16; i < 19; i++)
      m[i] = SIG1(m[i - 2]) + m[i - 7] + SIG0(m[i - 15]) + m[i - 16];

#undef ROT
#undef SHR

    clock_gettime(CLOCK_MONOTONIC, &t_start);
    {
      uint64_t chunks[3] = {
          0x4141410041414141,
          0x6161610061616161,
          0x2f2f2f002f2f2f2f,
      };

      for (uint8_t x = 0; x < chunks_n; x++) {
        uint8_t len = x == 2 ? 11 : 26;
        dim3 blocks(len, len, len);
        dim3 threads(len, len, 1);

        batch_size += len * len * len * len * len;

        uint64_t chunk = chunks[x];
        kernel<<<blocks, threads>>>(a, b, c, d, e, f, g, h, m[18],
                                    *(uint32_t *)(&msg[44]), chunk, block_cu,
                                    target_cu[x], mutex_cu[x]);

        cudaNoError(cudaPeekAtLastError());
      }

      for (uint8_t x = 0; x < chunks_n; x++) {
        uint32_t new_target;
        cudaNoError(cudaMemcpy(&new_target, target_cu[x], sizeof(uint32_t),
                               cudaMemcpyDeviceToHost));
        if (new_target < target) {
          target = new_target;
          cudaNoError(cudaMemcpy(block, block_cu, sizeof(uint8_t) * 12,
                                 cudaMemcpyDeviceToHost));
          printf("%08x ", target);
          for (uint8_t u = 0; u < 44; u += 4)
            printf("%c%c%c%c", msg[u + 3], msg[u + 2], msg[u + 1], msg[u + 0]);
          printf("%s\n", block);
        }
      }
    }

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    double time_taken = (t_end.tv_sec - t_start.tv_sec) +
                        (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;
    fprintf(stderr, "* %06x [%s] %3.2fs %3.2fGk/s | %08x => %s\n", batch,
            &msg[44], time_taken,
            batch_size / time_taken / 1000 / 1000 / 1000 * 26 * 26, target,
            block[0] > 0 ? (const char *)block : "<none>");
  }

  cudaNoError(cudaFree(block_cu));
  cudaNoError(cudaFree(target_cu));
  cudaNoError(cudaFree(mutex_cu));

  return 0;
}
